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
DP_OUTDIR_SUFFIX = "_DPoutput_1250iters_250burnin_seed123" # Change this between DPClust versions
#SNV_CLUSTERS_SUFFIX = "_1250iters_250burnin_bestClusterInfo.txt" #"_optimaInfo.txt" # Change this between DPClust versions
SNV_CLUSTERS_SUFFIX = "_optimaInfo.txt"
SNV_ASSIGNMENT_SUFFIX = "_1250iters_250burnin_bestConsensusAssignments.bed"
# Fraction of total SNVs assigned to a cluster to make it believable
FRAC_SNVS_CLUSTER = 0.01
# Minimum number of SNVs that make a cluster believable. This is jointly applied with the above filter, if either passed the cluster is kept
#MIN_NUM_SNVS_CLUSTER = 50 # Disabled for now as it may need to be increased, further investigation needed
CLONAL_PEAK_MIN_CCF = 0.9
CLONAL_PEAK_MAX_CCF = 1.1
PLOIDY_MAX_DIPLOID = 2.7
MEDIUM_SIZE_ABERRATION = 10 # CNA size in MB
LARGE_SIZE_ABERRATION = 50 # CNA size in MB
COVERAGE_FILES = c("../coverage/coverage_santa_cruz.txt", "../coverage/coverage_august_release_single.txt", "../coverage/coverage_august_release_multiple.txt", "../coverage/missing_coverage.tsv")
#GENDER_FILES = c("../gender/2015_05_15_santa_cruz_pilot_inferred_genders.txt", "../gender/2015_08_31_santa_cruz_pilot_inferred_genders_multiplesamples.txt", "../gender/2015_10_06_august_release_genders_single.txt", "../gender/2015_10_27_august_release_genders_multiple.txt", "../gender/2016_02_01_gender_unreleased_samples.txt", "../gender/2016_03_06_october_release_genders.txt")
GENDER_FILES = c("../gender/genders.txt")
ICGC_HISTOLOGY_FILE = "/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/clinical_data/pcawg_specimen_histology_August2016_v6_fixed.tsv"
ICGC_CLINICAL_FILE = "/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/clinical_data/pcawg_donor_clinical_August2016_v6_fixed.tsv"

#samplelist = read.table("2015_10_15_icgc_samples_pass.tsv", header=T, stringsAsFactors=F)
#samplelist = read.table("2015_10_29_icgc_samples_pass.tsv", header=T, stringsAsFactors=F)
sampleinput = read.table("samplenames_included.tsv", header=T, stringsAsFactors=F)
#is_refit = read.table("2015_10_29_sample_refitted_complete.tsv", header=T, stringsAsFactors=F)

#############################################################################################################################
# Annotate the cancer type
#############################################################################################################################
projectcode = sampleinput$projectcode
cancer_type = unlist(lapply(sampleinput$projectcode, function(x) { unlist(strsplit(x, "-"))[1] } ))
samplelist = sampleinput$sampleid
print("Basic table")
output = data.frame(projectcode=projectcode, cancer_type=cancer_type, samplename=samplelist)

#############################################################################################################################
# Purity
#############################################################################################################################
getPurity = function(samplename) {
  cancer_type = output[output$samplename==samplename, ]$projectcode
  filename = paste(PATH_TO_BB, cancer_type, "/", samplename, "/", samplename, RHO_PSI_SUFFIX, sep="")
  if (file.exists(filename)) {
  	  return(read.table(filename, header=T, stringsAsFactors=F)["FRAC_GENOME", "rho"])
  } else {
	  return(NA)
  }
}
print("Purity")
purity = unlist(lapply(samplelist, getPurity))
output$purity = purity

#############################################################################################################################
# Ploidy
#############################################################################################################################
getPloidy = function(samplename) {
  cancer_type = output[output$samplename==samplename, ]$projectcode
  subclones = read.table(paste(PATH_TO_BB, cancer_type, "/", samplename, "/", samplename, SUBCLONES_SUFFIX, sep=""), header=T, stringsAsFactors=F)
  subclones$length = round((subclones$endpos-subclones$startpos)/1000)
  cn_state_one = (subclones$nMaj1_A+subclones$nMin1_A)*subclones$frac1_A
  cn_state_two = ifelse(!is.na(subclones$frac2_A), (subclones$nMaj2_A+subclones$nMin2_A)*subclones$frac2_A, 0)
  ploidy = sum((cn_state_one+cn_state_two) * subclones$length, na.rm=T) / sum(subclones$length, na.rm=T)
  return(ploidy)
}

print("Ploidy")
ploidy = unlist(lapply(samplelist, getPloidy))
output$ploidy = round(ploidy, 3)

#############################################################################################################################
# num and frac clonal and num_subclones
#############################################################################################################################
getSubclonesAndAssignments = function(samplename, min_clonal_ccf=CLONAL_PEAK_MIN_CCF, max_clonal_ccf=CLONAL_PEAK_MAX_CCF) {
  sample_dp_dir = paste(PATH_TO_DP, samplename, DP_OUTDIR_SUFFIX, "/", sep="")
  # Make sure the files exist, this does seem to happen every once in a while
  if (!file.exists(paste(sample_dp_dir, samplename, SNV_CLUSTERS_SUFFIX, sep="")) | !file.exists(paste(sample_dp_dir, samplename, SNV_ASSIGNMENT_SUFFIX, sep=""))) {
    if (!file.exists(paste(sample_dp_dir, samplename, SNV_CLUSTERS_SUFFIX, sep=""))) {
      warning(paste(samplename, "no clusters"))
    }
    
    if (!file.exists(paste(sample_dp_dir, samplename, SNV_ASSIGNMENT_SUFFIX, sep=""))) {
      warning(paste(samplename, "no assignments"))
    }
    return(list(NA, NA, NA, NA, NA))
  }
  
  clusters = read.table(paste(sample_dp_dir, samplename, SNV_CLUSTERS_SUFFIX, sep=""), header=T, stringsAsFactors=F)
  assignments = table(read.table(paste(sample_dp_dir, samplename, SNV_ASSIGNMENT_SUFFIX, sep=""), header=T, stringsAsFactors=F)$cluster)
  total_muts = sum(assignments)
#  kept_clusters = names(assignments)[assignments > (total_muts*FRAC_SNVS_CLUSTER) | assignments > MIN_NUM_SNVS_CLUSTER]
#  kept_clusters = names(assignments)[assignments > (total_muts*FRAC_SNVS_CLUSTER)]
  kept_clusters = names(assignments)
  
  # Count the number of subclones and SNVs assigned
  num_clonal = 0
  num_subclonal = 0
  num_subclones = 0
  num_superclones = 0
  num_superclonal = 0
  cluster_locations = c()
  cluster_sizes = c()
  for (cluster in rev(kept_clusters)) {
    if (clusters[clusters$cluster.no==cluster,]$location > min_clonal_ccf) {
	   # Clonal
	   num_clonal = num_clonal + assignments[cluster]

	   if (clusters[clusters$cluster.no==cluster,]$location > max_clonal_ccf) {
		   # Superclonal
		   num_superclones = num_superclones + 1
		   num_superclonal = num_superclonal + assignments[cluster]
	   }
    } else {
      # Subclonal
      num_subclonal = num_subclonal + assignments[cluster]
      num_subclones = num_subclones + 1
    }
    cluster_locations = c(cluster_locations, clusters[clusters$cluster.no==cluster,]$location)
    cluster_sizes = c(cluster_sizes, assignments[cluster])
  }
  return(list(clust_stats=list(num_subclones, num_clonal, num_subclonal, num_superclones, num_superclonal), clust_details=list(cluster_locations, cluster_sizes)))
}
print("Num subclones")
res = lapply(samplelist, getSubclonesAndAssignments)
clust_stats = lapply(res, function(x) x$clust_stats)
# Saving these for later in the script to be appended to the table
clust_details = lapply(res, function(x) x$clust_details)

res = as.data.frame(matrix(unlist(clust_stats), ncol=5, byrow=T))
colnames(res) = c("num_subclones", "num_clonal", "num_subclonal", "num_superclones", "num_superclonal")
res$frac_clonal = round(res$num_clonal / (res$num_subclonal+res$num_clonal), 3)
output = data.frame(output, res)

#############################################################################################################################
# noCNA clonal  subclonal 
#############################################################################################################################

# normalCN means what should be considered the normal CN status of each allele
# Allows for discrimination between diploid and tetraploid samples
getCNstatus = function(segment, normalCN=1) {
  if (is.na(segment$nMaj1_A) | is.na(segment$nMaj1_A)) {
     status = "noCNA"
  } else if (segment$frac1_A==1) {
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

getCNAFractions = function(samplename, data_table, max_ploid_diploid=PLOIDY_MAX_DIPLOID, medium_size_aberration=MEDIUM_SIZE_ABERRATION, large_size_aberration=LARGE_SIZE_ABERRATION) {
  ploidy = data_table[data_table$samplename==samplename,]$ploidy
  cancer_type = data_table[data_table$samplename==samplename,]$projectcode
  subclones = read.table(paste(PATH_TO_BB, cancer_type, "/", samplename, "/", samplename, SUBCLONES_SUFFIX, sep=""), header=T, stringsAsFactors=F)
  subclones$len = round((subclones$endpos-subclones$startpos)/1000000)
  genome_len = sum(as.numeric(subclones$len))

  # Use ploidy to determine if there has been a WGD
  # If there is we assume 2 is the normal CN state and reason from there
  cn_sample = list(noCNA=0, subclonal=0, clonal=0)
  ploid = ifelse(ploidy<max_ploid_diploid, "diploid", "tetraploid")
  # Check if diploid or tetraploid
  if (ploidy <= max_ploid_diploid) {
    # diploid
    normalCN = 1
  } else {
    # tetraploid
    normalCN = 2
  }
  
  chroms = unique(subclones$chr)
  for (chrom in chroms) {
    subclones.c = subclones[subclones$chr==chrom,]
    for (i in 1:nrow(subclones.c)) {
      status = getCNstatus(subclones.c[i,], normalCN)

      # Add the fraction of the genome of this segment to the total of the category determined right above (status variable)
      frac = subclones.c[i,]$len
      cn_sample[names(cn_sample)==status][[1]] = cn_sample[names(cn_sample)==status][[1]] + frac
    }
  }
  
  has_clonal_cna = F
  has_clonal_cna_medium = F
  has_clonal_cna_large = F
  if (any(subclones$frac1_A==1 & !(subclones$nMaj1_A==normalCN & subclones$nMaj1_A==subclones$nMin1_A))) {
    has_clonal_cna = T
  }
  
  if (any(subclones$frac1_A==1 & !(subclones$nMaj1_A==normalCN & subclones$nMaj1_A==subclones$nMin1_A) & subclones$len > medium_size_aberration)) {
    has_clonal_cna_medium = T
  }
  
  if (any(subclones$frac1_A==1 & !(subclones$nMaj1_A==normalCN & subclones$nMaj1_A==subclones$nMin1_A) & subclones$len > large_size_aberration)) {
    has_clonal_cna_large = T
  }

  return(list(noCNA=round(cn_sample["noCNA"][[1]]/genome_len, 3), 
              clonal=round(cn_sample["clonal"][[1]]/genome_len, 3), 
              subclonal=round(cn_sample["subclonal"][[1]]/genome_len, 3), 
              ploidy_category=ploid,
              has_clonal_cna=has_clonal_cna,
              has_clonal_cna_medium=has_clonal_cna_medium,
              has_clonal_cna_large=has_clonal_cna_large))
}
print("Copy number")
res = as.data.frame(matrix(unlist(
  lapply(samplelist, function(x, ploidy) { getCNAFractions(x, output) })), 
  ncol=7, byrow=T))
colnames(res) = c("noCNA", "clonal", "subclonal", "ploidy_category", "has_clonal_cna", "has_clonal_cna_medium", "has_clonal_cna_large")
# Only add the required columns here
output = data.frame(output, res[, !grepl("has_clonal", colnames(res))])

# Save these for adding later
clonal_cna_inventory = res[, grepl("has_clonal", colnames(res))]

#############################################################################################################################
# coverage 
#############################################################################################################################
coverage = NULL
for (infile in COVERAGE_FILES) {
  d = read.table(infile, header=F, stringsAsFactors=F)
  colnames(d) = c("tumour", "cov_tumour", "normal", "cov_normal")
  coverage = rbind(coverage, d)
}
print("Coverage")
row_match = match(output$samplename, coverage$tumour)
output$cov_tumour = coverage[row_match,]$cov_tumour
output$cov_normal = coverage[row_match,]$cov_normal

#############################################################################################################################
# Sex
#############################################################################################################################
sex = NULL
for (infile in GENDER_FILES) {
  d = read.table(infile, header=T, stringsAsFactors=F)
  sex = rbind(sex, d)
}
print("Sex")
row_match = match(output$samplename, sex$tumour)
output$inferred_sex = sex[row_match,]$pred_gender

#############################################################################################################################
# power calculation
#############################################################################################################################
print("Power")
output$nrpcc = round((output$purity) / (output$purity*output$ploidy + (1-output$purity)*2) * output$cov_tumour, 3)

# Original implementation
#output$nrpcc3 = output$cov_tumour/output$ploidy*output$purity

#############################################################################################################################
# Add column whether a sample has been refit
#############################################################################################################################
#row_match = match(output$samplename, is_refit$samplename)
#output$refit = is_refit[row_match,2]
output$refit = NA

#############################################################################################################################
# Annotations from ICGC, ids and histology
#############################################################################################################################
print("ICGC Annotations")
anno = readr::read_tsv(ICGC_HISTOLOGY_FILE)
annotations = lapply(1:nrow(sampleinput), function(i) {
  sel = sampleinput$icgc_sample_id[i]==anno$icgc_sample_id & sampleinput$icgc_donor_id[i]==anno$icgc_donor_id & anno$specimen_library_strategy=="WGS"
  if (sum(sel) > 0) {
    anno[sel, c("histology_abbreviation", "icgc_sample_id", "icgc_donor_id")]
  } else {
    data.frame(histology_abbreviation=NA, icgc_sample_id=NA, icgc_donor_id=NA)
  }})
annotations = do.call(rbind, annotations)
output = data.frame(output, annotations)

anno = readr::read_tsv(ICGC_CLINICAL_FILE)
annotations = lapply(1:nrow(sampleinput), function(i) {
  sel = sampleinput$icgc_donor_id[i]==anno$icgc_donor_id
  if (sum(sel) > 0) {
    anno[sel, c("donor_sex")]
  } else {
    data.frame(donor_sex=NA)
  }})
annotations = do.call(rbind, annotations)
colnames(annotations) = c("reported_sex")
output = data.frame(output, annotations)

#############################################################################################################################
# Add cluster details, unraveled earlier when collecting the cluster stats
#############################################################################################################################
print("Cluster details")
res = lapply(clust_details, function(x) data.frame(paste(round(x[[1]], 4), collapse=";"), paste(x[[2]], collapse=";")))
res = do.call(rbind, res)
colnames(res) = c("cluster_locations", "cluster_sizes")
output = data.frame(output, res)

#############################################################################################################################
# Add clonal CNA inventory
#############################################################################################################################
print("Clonal CNA inventory")
output = data.frame(output, clonal_cna_inventory)

#############################################################################################################################
# save output
#############################################################################################################################
write.table(output, file="summary_table_2.txt", sep="\t", row.names=F, quote=F)

warnings()
q(save="no")
