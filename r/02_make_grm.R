# =============================================================================
# 02_make_grm.R
#
# Create the genomic relatedness matrix (GRM) for GBLUP.
#
# Input:
#   - Output from prepare_data()
#
# Output:
#   - GRM
#   - inverse GRM
#   - paths to GRM files
#
# This file is specific to GRM-based models such as GBLUP.
# =============================================================================


make_analysis_grm <- function(
    prepared,
    ncores = 8,
    mem = 8 * 6000
) {

  # ---------------------------------------------------------------------------
  # 1. Extract prepared objects
  # ---------------------------------------------------------------------------

  pheno_data <- prepared$pheno_data
  geno_files <- prepared$geno_files
  qc_filters <- prepared$qc_filters
  response <- prepared$response
  sys_name <- prepared$sys_name

  analysis_inds <- unique(pheno_data$ringnr)


  # ---------------------------------------------------------------------------
  # 2. Create raw GRM files
  # ---------------------------------------------------------------------------

  grm_files <- make_raw_grm(
    analysis_inds = analysis_inds,

    bfile = gsub(
      ".{4}$",
      "",
      geno_files[2]
    ),

    frq_file = geno_files[5],

    genorate_ind = qc_filters$genorate_ind,
    genorate_snp = qc_filters$genorate_snp,

    ncores = ncores,
    mem = mem,

    maf = qc_filters$maf,

    response = response,
    geno_set = paste0(sys_name, "_70K")
  )


  # ---------------------------------------------------------------------------
  # 3. Load GRM and compute inverse GRM
  # ---------------------------------------------------------------------------

  grm_obj <- compute_grm_obj(
    frq_file = geno_files[5],
    rel_file = grm_files[2],
    id_file = grm_files[1],
    bim_file = grm_files[3],

    pheno_data = pheno_data,

    id_col = 1
  )


  # ---------------------------------------------------------------------------
  # 4. Basic checks / summary
  # ---------------------------------------------------------------------------

  message("GRM construction complete.")
  message(
    "GRM dimensions: ",
    nrow(grm_obj$grm),
    " x ",
    ncol(grm_obj$grm)
  )

  message(
    "Inverse GRM dimensions: ",
    nrow(grm_obj$inv_grm),
    " x ",
    ncol(grm_obj$inv_grm)
  )


  # ---------------------------------------------------------------------------
  # 5. Return GRM objects
  # ---------------------------------------------------------------------------

  list(
    grm = grm_obj$grm,
    inv_grm = grm_obj$inv_grm,
    add_val = grm_obj$add_val,
    grm_files = grm_files
  )
}
