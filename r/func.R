# =============================================================================
# Utility functions
# =============================================================================

# gets iid from .fam
get_genotyped_inds <- function(fam_file,
                               sel = 2 # 1: fid; 2: iid
)  {
  fread(fam_file, select = sel, data.table = FALSE, header = FALSE)[, 1]
}

do_qc <- function(fam_file,
                  ncores,
                  mem,
                  qc_filt,
                  keep_inds,
                  sys,
                  resp) {
  # Extract filepath without ".fam"
  file_root <- gsub(pattern = ".fam", replacement = "", x = fam_file)
  # Specify directory where files will be stored
  dir <- paste0("data/qc", "_", resp, "_", sys)
  # Load fam file
  fam <- fread(fam_file, select = c(1, 2), data.table = FALSE, header = FALSE)
  # Specify which individuals in fam file we want to keep
  fam_keep <- fam[fam$V1 %in% keep_inds, ]
  # Exclude samples with high heterozygosity
  fam_keep <- fam_keep[!grepl(pattern = ".*HIGHHET.*", x = fam_keep$V2), ]
  # Exclude samples with mismatches in sex
  fam_keep <- fam_keep[!grepl(pattern = ".*MISSEX.*", x = fam_keep$V2), ]

  # For inds. genotyped multiple times, keep genotyping
  fam_keep <- fam_keep[!duplicated(fam_keep$V1, fromLast = TRUE), ]

  # Create folder where files will be saved
  dir.create(dir, showWarnings = FALSE)

  # Create new fam file with only inds. we want to keep
  write.table(fam_keep,
              file = paste0(dir, "/keep.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  # Call plink
  exit_code <-
    system2(get_plink_path(),
            paste0("--bfile ", file_root, " ",
                   "--make-bed ",
                   # Compute allele frequencies:
                   "--freq ",
                   # Filter by minor allele frequency:
                   "--maf ", qc_filt$maf, " ",
                   # Filter SNPs by call rate:
                   "--geno ", qc_filt$genorate_snp, " ",
                   # Filter inds by call rate:
                   "--mind ", qc_filt$genorate_ind, " ",
                   "--chr-set 32 ", # Sparrow chromosomes
                   "--memory ", mem, " ",
                   "--keep ", dir, "/keep.txt ",
                   "--threads ", ncores, " ",
                   "--out ", dir, "/qc"))

  # Throw error if unsuccessful
  if (exit_code != 0) {
    stop("Error in plink")
  }

  # Return path to new files if successful
  paste0(dir, "/", c("keep.txt", paste0("qc.", c("fam", "bim", "bed", "frq"))))
}

pheno_wrangle <- function(filepath,
                          genotyped_inds,
                          islands,
                          y_col_name,
                          testing = NULL) {
  # Load data
  dat <- fread(file = filepath, header = TRUE, data.table = FALSE) %>%
    # Keep only the...
    dplyr::filter(ringnr %in% genotyped_inds) %>% # genotyped individuals,
    dplyr::filter(!is.na(get(y_col_name))) %>% # phenotyped individuals,
    dplyr::filter(locality %in% islands) # and measurements in the system
  # Merge Lurøy and Onøy
  dat$first_locality <- ifelse(dat$first_locality %in% c(331, 332),
                               33,
                               dat$first_locality)

  # Do some checks
  stopifnot(all(dat$adult_sex %in% c(1, 2)))

  # Create phenotype data frame, create factors where relevant
  dat <- data.frame(ringnr = dat$ringnr,
                    sex = factor(ifelse(dat$adult_sex == 1, "m", "f")),
                    y = dat[, y_col_name],
                    year = factor(dat$year),
                    month = factor(dat$month),
                    island = factor(dat$locality),
                    hatch_year = factor(dat$hatch_year),
                    first_island = factor(dat$first_locality),
                    age = dat$year - dat$hatch_year,
                    day_session = interaction(dat$ringnr,
                                              dat$year,
                                              dat$month,
                                              dat$day,
                                              drop = TRUE),
                    id1 = match(dat$ringnr, unique(dat$ringnr)),
                    id2 = match(dat$ringnr, unique(dat$ringnr)))

  # If testing code, subset to only 100 observations
  if (!is.null(testing)) {
    s <- sort(sample(dim(dat)[1], testing))
    while (any(summary(dat$island[s]) == 0)) {
      s <- sample(dim(dat)[1], testing)
    }
    dat <- dat[s, ]
    dat$id1 <- dat$id2 <- match(dat$ringnr, unique(dat$ringnr))
  }

  # Return data frame
  dat
}

make_raw_grm <- function(analysis_inds,
                         bfile,
                         frq_file,
                         test_islands = NULL,
                         train_islands = NULL,
                         ncores,
                         mem,
                         genorate_ind,
                         genorate_snp,
                         maf,
                         response,
                         geno_set,
                         rel_cutoff = 1 - 1e-8) {

  dir <- paste0("data/grm_", response, "_", geno_set)

  if (!is.null(train_islands)) {
    dir <- paste0(dir, "_train", train_islands$code)
  }
  if (!is.null(test_islands)) {
    dir <- paste0(dir, "_test", test_islands$code)
  }

  dir.create(dir, showWarnings = FALSE)

  fam <- fread(file = paste0(bfile, ".fam"), header = FALSE)

  write.table(fam[fam$V1 %in% analysis_inds, ],
              file = paste0(dir, "/keep.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  exit_code <- system2(get_plink_path(),
                       paste0("--bfile ", bfile, " ",
                              "--maf ", maf, " ",
                              "--keep ", dir, "/keep.txt ",
                              "--geno ", genorate_snp, " ",
                              "--mind ", genorate_ind, " ",
                              "--chr-set 32 ",
                              "--memory ", mem, " ",
                              "--read-freq ", frq_file, " ",
                              "--rel-cutoff ", rel_cutoff, " ",
                              "--make-rel square bin cov ", # calculate raw GRM
                              "--make-just-bim ", # to know which SNPS were kept
                              "--threads ", ncores, " ",
                              "--out ", dir, "/grm"))

  if (exit_code != 0) {
    stop("Error in plink")
  }

  paste0(dir, "/grm.", c("rel.id", "rel.bin", "bim"))
}

compute_grm_obj <- function(frq_file,
                            rel_file,
                            id_file,
                            bim_file,
                            pheno_data,
                            id_col) {

  snps <- fread(file = bim_file, select = 2, data.table = FALSE, header = FALSE)
  frq <- fread(file = frq_file, select = c("SNP", "MAF"), header = TRUE)
  frq_inc <- frq[frq$SNP %in% snps$V2]
  ids <- fread(file = id_file,
               select = id_col,
               data.table = FALSE,
               header = FALSE)[, 1]
  n <- length(ids)

  vr_grm <- matrix(readBin(rel_file,
                           "numeric",
                           n ^ 2),
                   nrow = n) *
    (dim(frq_inc)[1] - 1) / # convert from sample covariance
    (2 * sum(frq_inc$MAF * (1 - frq_inc$MAF))) # vanRaden scaling

  # Rename columns for INLA
  dimnames(vr_grm)[[1]] <- dimnames(vr_grm)[[2]] <-
    pheno_data[match(ids, pheno_data$ringnr), "id1"]

  # Check symmetry
  if (!isSymmetric(vr_grm))
    stop("GRM not symmetric")

  # Compute eigenvalues
  e_vals <- eigen(vr_grm, only.values = TRUE)$values
  # Add small value to diagonal to get pos. def. matrix
  add_val <- ifelse(tail(e_vals, 1) < 0, - tail(e_vals, 1), 0) + 1e-9
  vr_grm <- vr_grm + diag(add_val, nrow(vr_grm))

  # Check positive defitiveness:
  e_vals_new <- eigen(vr_grm, only.values = TRUE)$values
  if (!all(Im(e_vals_new) == 0 & Re(e_vals_new) > 0))
    stop("GRM not positive definite")

  lst(grm = vr_grm, inv_grm = solve(vr_grm), add_val)
}

make_prior <- function(pc_prec_upper_var,
                       var_init,
                       tau = 0.05,
                       ...) {

  # PC priors for random effects
  hyperpar_var <- list(
    prec = list(initial = log(1 / var_init),
                prior = "pc.prec",
                param = c(sqrt(pc_prec_upper_var), tau),
                fixed = FALSE))

  lst(hyperpar_var)
}

make_cv_test_sets <- function(analysis_inds,
                              num_folds = 10) {
  n <- length(analysis_inds)

  # Split individuals into n folds
  scrambled <- sample(n)
  folds <- cut(seq_len(n), num_folds)
  levels(folds) <- seq_len(num_folds)

  # Create list of test individuals in each fold
  lapply(1:num_folds, function(fo) {
    x <- analysis_inds[sort(scrambled[folds == fo])]
    attributes(x) <- list(fo = fo)
    x
  })
}

run_gp <- function(pheno_data,
                   train_inds,
                   test_inds, # The rest in analysis_inds are used for training
                   inverse_relatedness_matrix = NULL,
                   effects_vec,
                   prior,
                   verbose = TRUE,
                   control.compute.config = TRUE
) {

  n_test <- length(test_inds)
  n_train <- length(train_inds)
  n <- n_test + n_train

  inla_formula <- stats::reformulate(effects_vec, response = "y_na")

  pheno_data$y_na <- pheno_data$y
  pheno_data$y_na[pheno_data$ringnr %in% test_inds] <- NA

  INLA::inla(inla_formula,
             family = "gaussian",
             data = pheno_data,
             verbose = verbose,
             control.compute = list(config = control.compute.config),
             control.family = list(hyper = prior$hyperpar_var)) %>%
    INLA::inla.rerun() %>%
    INLA::inla.rerun()
}

inla_posterior_variances <- function(prec_marginal) {
  sigma_marg <- INLA::inla.tmarginal(function(x) 1 / x, prec_marginal)
  INLA::inla.zmarginal(sigma_marg, silent = TRUE)
}

extract_stats <- function(post_sample) {
  mean <- mean(post_sample)
  sd <- sd(post_sample)
  var <- var(post_sample)
  mode <- suppressWarnings(MCMCglmm::posterior.mode(post_sample))
  median <- median(post_sample)
  hpd_lower <- hdi(post_sample)["lower"]
  hpd_upper <- hdi(post_sample)["upper"]
  lst(mean, mode, median, sd, var, hpd_lower, hpd_upper)
}


run_gp_time <- function(pheno_data,
                        train_rows,
                        test_rows,
                        inverse_relatedness_matrix = NULL,
                        effects_vec,
                        prior,
                        verbose = TRUE,
                        control.compute.config = TRUE
) {

  keep_rows <- train_rows | test_rows
  new_data <- pheno_data[keep_rows, ]

  new_data$y_na <- new_data$y
  new_data$y_na[test_rows[keep_rows]] <- NA

  inla_formula <- stats::reformulate(effects_vec, response = "y_na")

  model <- INLA::inla(inla_formula,
             family = "gaussian",
             data = new_data,
             verbose = verbose,
             control.compute = list(config = control.compute.config),
             control.family = list(hyper = prior$hyperpar_var)) %>%
    INLA::inla.rerun() %>%
    INLA::inla.rerun()

  list(
    model = model,
    data = new_data,
    test_rows = which(new_data$year_num == unique(pheno_data$year_num[test_rows]))
  )
}
