# =============================================================================
# Introductory script with basics for house sparrow quant. gen.
# Author: Kenneth Aase
# Expand on or change the code for your own workflow/analysis,
# this is just meant to be a simple and instructive starting point.
# ====================================s=========================================
# =============================================================================
#                            Part 0: setup
# =============================================================================
# ---- Packages necesarry to run this script ----
library(INLA) # Cannot be installed from CRAN, for install instructions for your
# system see: https://www.r-inla.org/download-install
library(dplyr)
library(data.table)
library(MCMCglmm)
library(HDInterval)
library(Rfast)
library(tictoc)
# ---- Data paths ----
# File path to morphological phenotype data
pheno_file <- "data/AdultMorphology_20241117.csv"
# Check if file exist
stopifnot(file.exists(pheno_file))

# File paths to genomic data
orig_geno_files <-
  paste0("data/combined_200k_70k_sparrow_genotype_data/",
         "combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05",
         c(".map", ".ped", ".fam", ".bim", ".bed"))
# Check if files exist
for (file in orig_geno_files) {
  stopifnot(file.exists(file))
}
# ---- File path to plink program here ----
# The script relies on external calls to the very useful PLINK 1.9 program for
# genomic data handling. Download the version matching your system from here:
# https://www.cog-genomics.org/plink/1.9/
# And make a path to the program file:
get_plink_path <- function() {
  "/home/stefamu/PLINK/plink"
}
# Check if program file exists
stopifnot(file.exists(get_plink_path()))

# ----  Load some utility functions ----
# Some details are deliberately hidden away inside these functions, so if you
# want to make changes somewhere, it would be wise to check inside these
source("code/func.R")
set.seed(2) # Set seed for reproducibility
# =============================================================================
#                     Part 1: Quality Control (QC)
# =============================================================================

# Vector of ID codes (ring numbers) for genotyped individuals
genotyped_inds <- get_genotyped_inds(fam_file = orig_geno_files[3], sel = 1)
# number of genotyped individuals
length(genotyped_inds)


# Basic quality control of the genomic data to select individuals and SNPs
# Remove SNPs and individuals with too much missing data,
# and SNPs with too little genetic variation (maf)
qc_filters <- list(genorate_ind = 0.05,
                   genorate_snp = 0.1,
                   maf = 0.01)
#                   "indep-pairwise" = c(50,5,0.2))
# See inside the do_qc function for more details.
# For more filtering options, see https://www.cog-genomics.org/plink/1.9/filter
qc_overall <- do_qc(fam_file = orig_geno_files[3],
                    ncores = 8,
                    mem = 8 * 6000,
                    qc_filt = qc_filters,
                    keep_inds = genotyped_inds,
                    sys = "",
                    resp = "overall")
genotyped_inds_qc <- get_genotyped_inds(fam_file = qc_overall[2], sel = 1)
# No of individuals that pass quality control
length(genotyped_inds_qc)

# =============================================================================
#     Part 2:  Use the quality controlled data and investigate LD patterns for the Helgeland system
# =============================================================================
# The question was whether LD patterns look different if only looking at the Helgeland system alone (it turns out, they don't really)

ncores=8

# where to store the results
dir <- paste0("data/qc_all_helgeland")

if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

# We use the pedigree to select individuals
isls <- c(20, 22, 23, 24, 26, 27, 28, 33, 331, 332, 34, 35, 38)
sys_name <- "helgeland"

# read the Helgeland pedigree and extract the ids of the individuals in the system
d.pedigree <-fread("data/pedigree/helgeland_inds_FID_IID.txt",data.table=FALSE)
genotyped_inds_helgeland <- as.character(d.pedigree[,1])

# Create new plink files that contain the individuals in the Helgeland pedigree and that passed quality control. All responses (no subsetting for whether phenotypes were available)
geno_files_helgeland <- do_qc(fam_file = qc_overall[2],
      ncores = 8,
      mem = 8 * 6000,
      qc_filt = qc_filters,
      keep_inds = genotyped_inds_helgeland,
      sys = sys_name,
      resp = "all")

genotyped_inds_helgeland <- get_genotyped_inds(fam_file = geno_files_helgeland[2], sel = 1)
# No of individuals that pass quality control
length(genotyped_inds_helgeland)

# LD checks: Calculate pairwise LD for all SNPs within a window of 1000 kb, this time only for the helgeland sparrows
file_root_helgeland <- "data/qc_all_helgeland/qc"

system2(get_plink_path(),
        paste0("--bfile ", file_root_helgeland, " ",
               #"--make-bed ",
               "--r2 ",
               "--ld-window-kb 1000 " ,
               "--ld-window 99999 ",
               "--ld-window-r2 0 " ,
               "--chr-set 32 ", # Sparrow chromosomes
               "--threads ", ncores, " ",
               "--out ", paste0(dir, "/qc_ld_helgeland_1000kb")
        )
)

ld <- read.table(paste0(dir, "/qc_ld_helgeland_1000kb.ld"), header=TRUE)

# Visualize the results
hist(ld$R2)
hist(ld$R2[ld$R2>0.1])
sort(ld$R2,decreasing=TRUE)[1:500]

ld$midpoint <- (ld$BP_A + ld$BP_B)/2
# Creating bins of a certain width (in base pairs)
ld$chr_window <- floor(ld$midpoint / 500000)

window_ld <- ld %>%
  group_by(CHR_A, chr_window) %>%
  summarise(mean_r2 = mean(R2))

plot(window_ld$mean_r2,type="l")


# =============================================================================
#     Part 3: Do LD pruning for the Helgeland system
# =============================================================================

# remove SNPs with r² > 0.8 within 1000 kb windows
system2(get_plink_path(),
        paste0("--bfile ", file_root_helgeland, " ",
               "--indep-pairwise 1000kb 5 0.8 ",
               "--chr-set 32 ",
               "--threads ", ncores, " ",
               "--out ", paste0(dir, "/qc_helgeland_ldprune")
        )
)

# Now we actually filter the dataset to only keep the SNPs that remain after pruning and create new fam/bin/bam files
system2(get_plink_path(),
        paste0("--bfile ", file_root_helgeland, " ",
               "--extract ", paste0(dir, "/qc_helgeland_ldprune.prune.in "),
               "--make-bed ",
               "--chr-set 32 ",
               "--threads ", ncores, " ",
               "--out ", paste0(dir, "/qc_helgeland_ldpruned_data")
        )
)

file_root_pruned <- paste0(dir,"/qc_helgeland_ldpruned_data")

system2(get_plink_path(),
        paste0("--bfile ", file_root_pruned, " ",
               #"--make-bed ",
               "--r2 ",
               "--ld-window-kb 1000 " ,
               "--ld-window 99999 ",
               "--ld-window-r2 0 " ,
               "--chr-set 32 ", # Sparrow chromosomes
               "--threads ", ncores, " ",
               "--out ", paste0(dir, "/qc_ld_helgeland_pruned_1000kb")
        )
)

ld <- read.table(paste0(dir, "/qc_ld_helgeland_pruned_1000kb.ld"), header=TRUE)

# Visualize the results
hist(ld$R2)
hist(ld$R2[ld$R2>0.1])
sort(ld$R2,decreasing=TRUE)[1:500]

ld$midpoint <- (ld$BP_A + ld$BP_B)/2
# Creating bins of a certain width (in base pairs)
ld$chr_window <- floor(ld$midpoint / 500000)

window_ld <- ld %>%
  group_by(CHR_A, chr_window) %>%
  summarise(mean_r2 = mean(R2))

plot(window_ld$mean_r2,type="l")


# =============================================================================
#     Part 4: Use the pruned Helgeland data and do the SVD
# =============================================================================
########################
# Part 5.1: Select individuals that should be retained in the PCR:
########################

##### Just body mass measurements
response_colname <- "body_mass"
response <- "mass"

# For quality check, look up how many individuals are in the phenotype file
pheno_all <- fread(file = pheno_file, header = TRUE, data.table = FALSE)
dim(pheno_all)
length(unique(pheno_all$ringnr))

# ---- Create data frame of chosen phenotype data ----
pheno_data <- pheno_wrangle(
  filepath = pheno_file,
  genotyped_inds = genotyped_inds_helgeland,
  islands = isls,
  y_col_name = response_colname,
  # Use this to include only e.g. 250 observations, for fast testing of code
  # (can cause numerical problems if set too low):
  # testing = 2000
  # Use this to include all observations:
  testing = NULL
)
# Inspect:
head(pheno_data)
dim(pheno_data)
# No of unique individuals:
length(unique(pheno_data$ringnr))

# Also add the
froh <- data.table::fread(file = "data/FROH2.5_helgeland.txt")
pheno_data$f_roh = froh$FROH2.5[match(pheno_data$ringnr, froh$FID)]

########################
### Part 5.2: Do PCA in plink
########################

# Vector of ringnrs to include. Example:
analysis_inds <- unique(pheno_data$ringnr)

max_num_pc <- length(analysis_inds)

mem <- 4000 # TODO: memory in mb
ncores <- 4 # TODO: CPUs to use
bfile <- file_root_Helgeland_pruned <- paste0(dir,"/qc_helgeland_ldpruned_data")

# Read the fam file in order to make some checks
genotyped.inds.tmp <- get_genotyped_inds(fam_file = paste0(bfile,".fam"), sel = 1)
length(genotyped.inds.tmp)
# How many of the individuals in the pheno file are in the genotype file? It should be all:
sum(analysis_inds %in% genotyped.inds.tmp)
# -> That seems correct.
fam <- read.table(paste0(bfile, ".fam"), stringsAsFactors = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")
head(fam)

# Write out the first two rows of the fam file if the ID matches with the individuals that we have to analyze:
keep <- fam[fam$FID %in% analysis_inds, c("FID", "IID")]

# Create directory to store PCA files
dir <- paste0(bfile, "_pca_mass") # Directory name
dir.create(dir, showWarnings = FALSE)
# Create file containing ringnr of inds to include
# Must equal the first 2 columns in bfile.fam
write.table(
  keep,
  file = file.path(dir, "keep.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

keep_file <- file.path(dir, "keep.txt")

# Doing the SVD with PLINK
exit_code <- system2(
  get_plink_path(),
  paste0("--bfile ", bfile, " ",
         "--pca ", max_num_pc, " \'header\' ", # Do PCA with a certain no. PCs
        # "--maf 0.01 ", # QC: Filter by minor allele frequency
        "--keep ", dir, "/keep.txt ", # Include only the inds. in keep.txt
        # "--geno 0.1 ", # QC: Filter SNPs by call rate
        # "--mind 0.05 ", # QC: Filter inds. by call rate
         "--chr-set 32 ", # Sparrow chromosomes
         "--memory ", mem, " ",
         "--threads ", ncores, " ",
         "--out ", dir, "/pca")) # PCA result saved here

if (exit_code != 0) {
  stop("Error in plink")
}

###### Make PC matrix
library(data.table)
# Load PCA results - PC matrix
pca <- fread(file = paste0(dir, "/pca.eigenvec"))




###### Variance analysis

# Load eigenvalues
eigenvals <- fread(file = paste0(dir, "/pca.eigenval"))
# Find amount of variance explained by each PC cumulatively
variance_proportion <- cumsum(eigenvals$V1) / sum(eigenvals$V1)
plot(variance_proportion)

# Find number of PCs needed to explain a given prop. of the genomic variation
min_var_explained <- 0.8 # Set to desired amount (but the no of PCs is anyway determined in other ways)
enough_variance_explained <- cumsum(variance_proportion > min_var_explained)
(num_pc <- which(enough_variance_explained == 1)) # Number of PCs required
var_explained_by_num_pc <- variance_proportion[num_pc]

######

# Be careful: Plink gives out standardized PCs! So we have to re-scale them to explain the proportion of variance given by the eigenvalue:

# Number of PCs
n_PCs <- 1000

# Create PC matrix -
zz <- data.frame(pca[, -c("FID", "IID")])[match(analysis_inds, pca$FID), ]
zz <- zz[,1:n_PCs]

# Scale the PCs back so that they explain variance proportional to the eigenvalues
z_scaled <- sweep(as.matrix(zz), 2, sqrt(eigenvals[[1]][1:n_PCs]), "*")

# Then divide by the SD of the first eigenvalue, so that PC1 has a variance of 1 etc.
z_scaled_final <- z_scaled / sqrt(var(z_scaled)[1])

# =============================================================================
#     Part 5: Do the BPCRR on a trait
# =============================================================================
#################
### Fixed prior (formula (4))
#################

# Summing the total variance of the PCs to use as fixed prior
PCvar <- colVars(z_scaled_final)
tot_PCvar <- sum(PCvar)

# As VA is set to 0.3 times the phenotypic variance to obtain the fixed prior
varA <- 5.2 * 0.3
u.prior.var <- varA/tot_PCvar

formula.sim  = y ~ sex + month + f_roh +
  f(hatch_year, model="iid",hyper=list(
    prec=list(initial=log(1), prior="pc.prec",param=c(1,0.1))
  )) +
  f(first_island, model="iid",hyper=list(
    prec=list(initial=log(1), prior="pc.prec",param=c(1,0.1))
  )) +
  f(id1, model = "z", Z = z_scaled_final,
                           hyper=list(
                             prec=list(initial=log(1/u.prior.var),
                                       fixed=TRUE
                                       # fixed=TRUE fixes the variance at u.prior.var
                                       # fixed=FALSE would give the default priors, but other priors can be specified as well.
                             )
                           )
)

# Now we are calling inla() - this call takes a bit of time (use tic() and toc() to measure):
tic()
model.sim = inla(formula=formula.sim, family="gaussian",
                 data=pheno_data,
                 control.family=list(hyper = list(theta = list(initial=log(1),
                                                               prior="pc.prec",
                                                               param=c(1,0.1)))),
                 control.compute=list(dic=F, config =TRUE,
                                      return.marginals=FALSE
                 ), # To be able to resample from the inla object, we need config=TRUE, but makes computation slower.
                 num.threads=10 # Set the number of cores for parallel computation. Use 1-2 less than what you have on your machine.
)
toc()

summary(model.sim)

#################################

Nanimals <-length(analysis_inds)
nsamples <- 1000
sample.sim <- inla.posterior.sample(n=nsamples,model.sim)

# Extract samples of the breeding values; in each entry of the list the breeding values of all animals in one sample are stored
samples.a <-  lapply(1:nsamples, function(ii) {
  sample.sim[[ii]]$latent[substring(rownames(sample.sim[[1]]$latent),1,3)=="id1"][1:Nanimals]
})

VAs.sim <- unlist(lapply(samples.a, var))

# Estimated VA
mean(VAs.sim)

# Posterior distribution of Va
hist(VAs.sim)
