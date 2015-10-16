# Colnames to assemble: 
# cancer_type     
# samplename      
# num_subclones   
# purity  
# ploidy  
# num_clonal      
# num_subclonal   
# frac_clonal     
# noCNA   
# clonal  
# subclonal       
# cov_tumour      
# cov_normal      
# ploidy_category 
# nrpcc

PATH_TO_BB = "../battenberg/"
PATH_TO_DP = "../dirichlet/"
SUBCLONES_SUFFIX = "_subclones.txt"
RHO_PSI_SUFFIX = "_rho_and_psi.txt"
DP_OUTDIR_SUFFIX = "_DPoutput_1250iters_250burnin"
SNV_CLUSTERS_SUFFIX = "_optimaInfo.txt"
SNV_ASSIGNMENT_SUFFIX = "_1250iters_250burnin_bestConsensusAssignments.bed"
# Fraction of total SNVs assigned to a cluster to make it believable
FRAC_SNVS_CLUSTER = 0.01
COVERAGE_FILES = c("../coverage/coverage_santa_cruz.txt", "../coverage/coverage_august_release_single.txt")

samplelist = read.table("2015_10_15_icgc_samples_pass.tsv", header=T, stringsAsFactors=F)

#############################################################################################################################
# Annotate the cancer type
#############################################################################################################################
cancer_type = unlist(lapply(samplelist$projectcode, function(x) { unlist(strsplit(x, "-"))[1] } ))
samplelist = samplelist$sampleid

output = data.frame(cancer_type=cancer_type, samplename=samplelist)

#############################################################################################################################
# Purity
#############################################################################################################################
getPurity = function(samplename) {
  return(read.table(paste(PATH_TO_BB, samplename, RHO_PSI_SUFFIX, sep=""), header=T, stringsAsFactors=F)["FRAC_GENOME", "rho"])
}
purity = unlist(lapply(samplelist, getPurity))
output$purity = purity

#############################################################################################################################
# Ploidy
#############################################################################################################################
getPloidy = function(samplename) {
  subclones = read.table(paste(PATH_TO_BB, samplename, SUBCLONES_SUFFIX, sep=""), header=T, stringsAsFactors=F)
  subclones$length = round((subclones$endpos-subclones$startpos)/1000)
  cn_state_one = (subclones$nMaj1_A+subclones$nMin1_A)*subclones$frac1_A
  cn_state_two = ifelse(!is.na(subclones$frac2_A), (subclones$nMaj2_A+subclones$nMin2_A)*subclones$frac2_A, 0)
  ploidy = sum((cn_state_one+cn_state_two) * subclones$length) / sum(subclones$length)
  return(ploidy)
}

ploidy = unlist(lapply(samplelist, getPloidy))
output$ploidy = ploidy

#############################################################################################################################
# num and frac clonal and num_subclones
#############################################################################################################################
getSubclonesAndAssignments = function(samplename, min_clonal_ccf=0.9) {
  sample_dp_dir = paste(PATH_TO_DP, samplename, DP_OUTDIR_SUFFIX, "/", sep="")
  # Make sure the files exist, this does seem to happen every once in a while
  if (!file.exists(paste(sample_dp_dir, samplename, SNV_CLUSTERS_SUFFIX, sep="")) | !file.exists(paste(sample_dp_dir, samplename, SNV_ASSIGNMENT_SUFFIX, sep=""))) {
    return(list(NA, NA, NA))
  }
  
  clusters = read.table(paste(sample_dp_dir, samplename, SNV_CLUSTERS_SUFFIX, sep=""), header=T, stringsAsFactors=F)
  assignments = table(read.table(paste(sample_dp_dir, samplename, SNV_ASSIGNMENT_SUFFIX, sep=""), header=T, stringsAsFactors=F)$cluster)
  total_muts = sum(assignments)
  kept_clusters = names(assignments)[assignments > total_muts*FRAC_SNVS_CLUSTER]
  
  # Count the number of subclones and SNVs assigned
  num_clonal = 0
  num_subclonal = 0
  num_subclones = 0
  for (cluster in kept_clusters) {
    if (clusters[clusters$cluster.no==cluster,]$location > min_clonal_ccf) {
      # Clonal
      num_clonal = num_clonal + assignments[cluster]
    } else {
      # Subclonal
      num_subclonal = num_subclonal + assignments[cluster]
      num_subclones = num_subclones + 1
    }
  }
  return(list(num_subclones, num_clonal, num_subclonal))
}

res = as.data.frame(matrix(unlist(lapply(samplelist, getSubclonesAndAssignments)), ncol=3, byrow=T))
colnames(res) = c("num_subclones", "num_clonal", "num_subclonal")
res$frac_clonal = round(res$num_subclonal / (res$num_subclonal+res$num_clonal), 3)
output = data.frame(output, res)

#############################################################################################################################
# noCNA clonal  subclonal 
#############################################################################################################################

# normalCN means what should be considered the normal CN status of each allele
# Allows for discrimination between diploid and tetraploid samples
getCNstatus = function(segment, normalCN=1) {
  if (segment$frac1_A==1) {
    # Check if normal or abberrated
    if (segment$nMaj1_A==normalCN & segment$nMaj1_A==segment$nMin1_A) {
      status = "noCNA"
    } else {
      status = "clonal"
    }
    
  } else {
    status = "subclonal"
  }
  return(status)
}

getCNAFractions = function(samplename, data_table, max_ploid_diploid=2.7) {
  ploidy = data_table[data_table$samplename==samplename,]$ploidy
  subclones = read.table(paste(PATH_TO_BB, samplename, SUBCLONES_SUFFIX, sep=""), header=T, stringsAsFactors=F)
  #dat = dat[,1:13]
  subclones$len = round((subclones$endpos-subclones$startpos)/1000)
  genome_len = sum(as.numeric(subclones$len))

  # Use ploidy to determine if there has been a WGD
  # If there is we assume 2 is the normal CN state and reason from there
  cn_sample = list(noCNA=0, subclonal=0, clonal=0)
  ploid = ifelse(ploidy<max_ploid_diploid, "diploid", "tetraploid")
  chroms = unique(subclones$chr)
  for (chrom in chroms) {
    subclones.c = subclones[subclones$chr==chrom,]
    for (i in 1:nrow(subclones.c)) {
      # Check if diploid or tetraploid
      if (ploidy <= max_ploid_diploid) {
        # diploid
        status = getCNstatus(subclones.c[i,], 1)
      } else {
        # tetraploid
        status = getCNstatus(subclones.c[i,], 2)
      }

      # Add the fraction of the genome of this segment to the total of the category determined right above (status variable)
      frac = subclones.c[i,]$len
      cn_sample[names(cn_sample)==status][[1]] = cn_sample[names(cn_sample)==status][[1]] + frac
    }
  }

  return(list(noCNA=round(cn_sample["noCNA"][[1]]/genome_len, 3), 
              clonal=round(cn_sample["clonal"][[1]]/genome_len, 3), 
              subclonal=round(cn_sample["subclonal"][[1]]/genome_len, 3), 
              ploidy_category=ploid))
}

res = as.data.frame(matrix(unlist(
  lapply(samplelist, function(x, ploidy) { getCNAFractions(x, output) })), 
  ncol=4, byrow=T))
colnames(res) = c("noCNA", "clonal", "subclonal", "ploidy_category")
output = data.frame(output, res)

#############################################################################################################################
# coverage 
#############################################################################################################################
coverage = NULL
for (infile in COVERAGE_FILES) {
  d = read.table(infile, header=F, stringsAsFactors=F)
  colnames(d) = c("tumour", "cov_tumour", "normal", "cov_normal")
  coverage = rbind(coverage, d)
}
row_match = match(output$samplename, coverage$tumour)
output$cov_tumour = coverage[row_match,]$cov_tumour
output$cov_normal = coverage[row_match,]$cov_normal

#############################################################################################################################
# power calculation
#############################################################################################################################
output$nrpcc = output$purity*output$cov_tumour/output$ploidy / (output$purity*output$cov_tumour/output$ploidy + (1-output$purity)*output$cov_tumour*2) * output$cov_tumour

#############################################################################################################################
# save output
#############################################################################################################################
write.table(output, file="summary_table.txt", sep="\t", row.names=F, quote=F)