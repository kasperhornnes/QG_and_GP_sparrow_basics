# =============================================================================
# 02_make_pcs.R
#
# Create the genomic principal-component representation used for BPCRR.
#
# This file:
#   1. LD-prunes the QC'ed genotype data
#   2. Performs PCA with PLINK
#   3. Rescales the PCs for use in BPCRR
#
# It does NOT fit the BPCRR model.
# =============================================================================


make_bpcrr_pcs <- function(
    prepared,
    n_pcs = 1000,
    ncores = 8,
    mem = 4000,
    out_dir = "data/bpcrr_pca"
) {

  pheno_data <- prepared$pheno_data
  geno_files <- prepared$geno_files

  # Order is important:
  # pheno_wrangle defines id1 according to unique(ringnr),
  # so we keep the same ordering here.
  analysis_inds <- unique(pheno_data$ringnr)


  # ---------------------------------------------------------------------------
  # 1. Create output directory
  # ---------------------------------------------------------------------------

  if (!dir.exists(out_dir)) {
    dir.create(
      out_dir,
      recursive = TRUE
    )
  }


  # Root of the final QC'ed PLINK files
  bfile <- sub(
    "\\.fam$",
    "",
    geno_files[2]
  )


  # ---------------------------------------------------------------------------
  # 2. LD pruning
  # ---------------------------------------------------------------------------

  prune_root <- file.path(
    out_dir,
    "ldprune"
  )

  exit_code <- system2(
    get_plink_path(),
    paste0(
      "--bfile ", bfile, " ",
      "--indep-pairwise 1000kb 5 0.8 ",
      "--chr-set 32 ",
      "--threads ", ncores, " ",
      "--memory ", mem, " ",
      "--out ", prune_root
    )
  )

  if (exit_code != 0) {
    stop("PLINK LD pruning failed.")
  }


  # ---------------------------------------------------------------------------
  # 3. Create LD-pruned genotype files
  # ---------------------------------------------------------------------------

  pruned_root <- file.path(
    out_dir,
    "ldpruned_data"
  )

  exit_code <- system2(
    get_plink_path(),
    paste0(
      "--bfile ", bfile, " ",
      "--extract ", prune_root, ".prune.in ",
      "--make-bed ",
      "--chr-set 32 ",
      "--threads ", ncores, " ",
      "--memory ", mem, " ",
      "--out ", pruned_root
    )
  )

  if (exit_code != 0) {
    stop("Creating LD-pruned genotype data failed.")
  }


  # ---------------------------------------------------------------------------
  # 4. Make keep file
  # ---------------------------------------------------------------------------

  fam <- read.table(
    paste0(pruned_root, ".fam"),
    stringsAsFactors = FALSE
  )

  colnames(fam)[1:2] <- c(
    "FID",
    "IID"
  )

  keep <- fam[
    as.character(fam$FID) %in%
      as.character(analysis_inds),
    c("FID", "IID")
  ]

  if (nrow(keep) != length(analysis_inds)) {
    stop(
      "Not all phenotype individuals were found ",
      "in the LD-pruned genotype data."
    )
  }

  keep_file <- file.path(
    out_dir,
    "keep.txt"
  )

  write.table(
    keep,
    file = keep_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 5. PCA with PLINK
  # ---------------------------------------------------------------------------

  # Cannot request more PCs than available individuals
  n_pcs <- min(
    n_pcs,
    length(analysis_inds) - 1
  )

  pca_root <- file.path(
    out_dir,
    "pca"
  )

  exit_code <- system2(
    get_plink_path(),
    paste0(
      "--bfile ", pruned_root, " ",
      "--pca ", n_pcs, " header ",
      "--keep ", keep_file, " ",
      "--chr-set 32 ",
      "--memory ", mem, " ",
      "--threads ", ncores, " ",
      "--out ", pca_root
    )
  )

  if (exit_code != 0) {
    stop("PLINK PCA failed.")
  }


  # ---------------------------------------------------------------------------
  # 6. Read PCA results
  # ---------------------------------------------------------------------------

  pca <- data.table::fread(
    paste0(pca_root, ".eigenvec"),
    data.table = FALSE
  )

  eigenvalues <- data.table::fread(
    paste0(pca_root, ".eigenval"),
    data.table = FALSE
  )

  # PLINK versions may call first column #FID
  names(pca)[1:2] <- c(
    "FID",
    "IID"
  )


  # ---------------------------------------------------------------------------
  # 7. Put PCs in same individual order as pheno_data / id1
  # ---------------------------------------------------------------------------

  idx <- match(
    as.character(analysis_inds),
    as.character(pca$FID)
  )

  if (anyNA(idx)) {
    stop(
      "Some analysis individuals are missing ",
      "from the PCA output."
    )
  }

  pc_raw <- as.matrix(
    pca[
      idx,
      -(1:2),
      drop = FALSE
    ]
  )


  # ---------------------------------------------------------------------------
  # 8. Rescale PCs
  #
  # PLINK outputs standardized PCs. Following the supplied BPCRR script,
  # multiply each PC by sqrt(eigenvalue), then scale relative to PC1.
  # ---------------------------------------------------------------------------

  z_scaled <- sweep(
    pc_raw,
    2,
    sqrt(eigenvalues[[1]][seq_len(n_pcs)]),
    "*"
  )

  z_scaled_final <-
    z_scaled / sqrt(stats::var(z_scaled)[1, 1])

  rownames(z_scaled_final) <-
    as.character(analysis_inds)


  # ---------------------------------------------------------------------------
  # 9. Summary
  # ---------------------------------------------------------------------------

  message("BPCRR PCA preparation complete.")
  message(
    "Individuals: ",
    nrow(z_scaled_final)
  )
  message(
    "PCs: ",
    ncol(z_scaled_final)
  )


  # ---------------------------------------------------------------------------
  # 10. Return
  # ---------------------------------------------------------------------------

  list(
    Z = z_scaled_final,
    eigenvalues = eigenvalues[[1]][seq_len(n_pcs)],
    analysis_inds = analysis_inds,
    n_pcs = n_pcs,
    pruned_bfile = pruned_root,
    pca_root = pca_root
  )
}
