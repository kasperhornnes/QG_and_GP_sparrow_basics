# =============================================================================
# 01_prepare_data.R
#
# Prepare phenotype and genotype data for genomic prediction analyses.
#
# This file:
#   1. Checks required files/programs
#   2. Performs overall genomic QC
#   3. Selects and cleans phenotype observations
#   4. Performs genomic QC for the final analysis individuals
#
# It does NOT:
#   - Create temporal train/test splits
#   - Fit prediction models
#   - Create the GRM
#   - Calculate prediction accuracy
# =============================================================================


prepare_data <- function(
    pheno_file,
    orig_geno_files,
    islands,
    sys_name,
    froh_file,
    response_colname = "body_mass",
    response = "mass",
    testing = NULL,
    qc_filters = list(
      genorate_ind = 0.05,
      genorate_snp = 0.10,
      maf = 0.01
    ),
    ncores = 8,
    mem = 8 * 6000
) {

  # ---------------------------------------------------------------------------
  # 1. Check input files
  # ---------------------------------------------------------------------------

  if (!file.exists(pheno_file)) {
    stop("Phenotype file not found: ", pheno_file)
  }

  missing_geno_files <- orig_geno_files[!file.exists(orig_geno_files)]

  if (length(missing_geno_files) > 0) {
    stop(
      "Missing genotype file(s):\n",
      paste(missing_geno_files, collapse = "\n")
    )
  }

  if (!file.exists(get_plink_path())) {
    stop("PLINK executable not found at: ", get_plink_path())
  }



  # ---------------------------------------------------------------------------
  # 2. Find all genotyped individuals
  # ---------------------------------------------------------------------------

  genotyped_inds <- get_genotyped_inds(
    fam_file = orig_geno_files[3],
    sel = 1
  )


  # ---------------------------------------------------------------------------
  # 3. Overall genotype QC
  # ---------------------------------------------------------------------------

  qc_overall <- do_qc(
    fam_file = orig_geno_files[3],
    ncores = ncores,
    mem = mem,
    qc_filt = qc_filters,
    keep_inds = genotyped_inds,
    sys = "",
    resp = "overall"
  )


  # Individuals remaining after QC
  genotyped_inds_qc <- get_genotyped_inds(
    fam_file = qc_overall[2],
    sel = 1
  )


  # ---------------------------------------------------------------------------
  # 4. Prepare phenotype data
  # ---------------------------------------------------------------------------

  pheno_data <- pheno_wrangle(
    filepath = pheno_file,
    genotyped_inds = genotyped_inds_qc,
    islands = islands,
    y_col_name = response_colname,
    testing = testing
  )


  # Numeric year is useful for all temporal analyses
  pheno_data$year_num <- as.numeric(
    as.character(pheno_data$year)
  )

  # Add genomic inbreeding coefficient
  froh <- data.table::fread(
    froh_file,
    data.table = FALSE
  )

  idx <- match(
    as.character(pheno_data$ringnr),
    as.character(froh$FID)
  )

  if (all(is.na(idx))) {
    stop("No phenotype IDs matched the FROH file.")
  }

  pheno_data$f_roh <- froh[["FROH2.5"]][idx]

  message(
    "Added f_roh: ",
    "f_roh" %in% names(pheno_data)
  )


  # ---------------------------------------------------------------------------
  # 5. Final genotype QC for individuals in this phenotype dataset
  # ---------------------------------------------------------------------------

  analysis_inds <- unique(pheno_data$ringnr)

  geno_files <- do_qc(
    fam_file = qc_overall[2],
    ncores = ncores,
    mem = mem,
    qc_filt = qc_filters,
    keep_inds = analysis_inds,
    sys = sys_name,
    resp = response
  )


  # ---------------------------------------------------------------------------
  # 6. Basic checks / summary
  # ---------------------------------------------------------------------------

  message("Data preparation complete.")
  message("Number of phenotype observations: ", nrow(pheno_data))
  message("Number of unique individuals: ", length(analysis_inds))
  message(
    "Years: ",
    min(pheno_data$year_num),
    " - ",
    max(pheno_data$year_num)
  )


  # ---------------------------------------------------------------------------
  # 7. Return common analysis data
  # ---------------------------------------------------------------------------

  list(
    pheno_data = pheno_data,
    analysis_inds = analysis_inds,
    geno_files = geno_files,
    qc_overall = qc_overall,
    qc_filters = qc_filters,
    response = response,
    response_colname = response_colname,
    sys_name = sys_name
  )
}
