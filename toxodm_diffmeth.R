#!/usr/bin/env Rscript
# ToxoDM phase-2 differential methylation -- Sarah's 4-model plan.
# Per comparison, four models (q-values reported, NO cutoff) + leave-one-out:
#   model_0      ~ status                      (unadjusted)
#   model_oligo  ~ status + Oligo_NNLS
#   model_micro  ~ status + Micro_NNLS
#   model_neuron ~ status + Neuron_total_NNLS
# Composition covariates read from the validated 4-type NNLS table.
# Self-contained: reads methylKit txt + proportions from disk, no session state.
#
# ART confound (per Sarah): A,B,E not on ART; C,D on ART.
#   comp1 = 2024neg(onART) vs 2026neg(offART)  -> characterizes ART/cohort, NOT toxo
#                                                  (both groups toxo-negative)
#   comp2 = 2024pos(offART) vs 2024neg(onART)  -> toxo contrast, ART-confounded
#   comp3 = 2024pos(offART) vs 2024neg(onART)  -> toxo contrast, ART-confounded (BG)
suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
})

MKDIR  <- "/projects1/ToxoDM/working/methylkit"
PROPF  <- "/projects1/ToxoDM/reports/nnls_proportions_4type.tsv"
OUTDIR <- "/projects1/ToxoDM/reports"
QCDIR  <- file.path(OUTDIR, "QC")
dir.create(QCDIR, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(OUTDIR, "diffmeth_run_comp4.log"), open = "wt")
sink(log_con, split = TRUE)      # split=TRUE: write to log AND console
sink(log_con, type = "message")  # also capture messages/warnings

prop <- fread(PROPF)
prop[, key := paste(cohort, tissue, sample, sep = "_")]

# the four models, as covariate-column-sets (NULL = unadjusted)
MODELS <- list(
  model_0      = NULL,
  model_oligo  = "Oligo_NNLS",
  model_micro  = "Micro_NNLS",
  model_neuron = "Neuron_total_NNLS"
)

# ---- comparison definitions (reference group = 0, treatment group = 1) ----------
NEG2024_FMC <- c("2024_FMC_Atpsy_15","2024_FMC_Atpsy_20","2024_FMC_Atpsy_21","2024_FMC_LG29")
NEG2024_BG  <- c("2024_BG_Atpsy_15","2024_BG_Atpsy_20","2024_BG_Atpsy_21","2024_BG_LG29")
POS_CC      <- c("2024_CC_CA192","2024_CC_CA346","2024_CC_CB239","2024_CC_CC249","2024_CC_CE167","2024_CC_CE234")
POS_BG      <- c("2024_BG_CA192","2024_BG_CA346","2024_BG_CB239","2024_BG_CC249","2024_BG_CE167","2024_BG_CE234")
NEG2026_CC  <- c("2026_CC_LG30","2026_CC_LG31","2026_CC_LG52")

comparisons <- list(
  comp1_ART_neg2024_vs_neg2026 = list(
    keys0 = NEG2024_FMC, keys1 = NEG2026_CC,
    label0 = "2024neg_onART_FMC", label1 = "2026neg_offART_CC",
    note   = "Both toxo-negative; contrast is ART/cohort/region, NOT toxo."
  ),
  comp2_toxo_pos_vs_neg_CC = list(
    keys0 = NEG2024_FMC, keys1 = POS_CC,
    label0 = "neg_onART_FMC", label1 = "pos_offART_CC",
    note   = "Toxo contrast, confounded with ART (pos off-ART, neg on-ART)."
  ),
  comp3_toxo_pos_vs_neg_BG = list(
    keys0 = NEG2024_BG, keys1 = POS_BG,
    label0 = "neg_onART_BG", label1 = "pos_offART_BG",
    note   = "Toxo contrast (BG), confounded with ART (pos off-ART, neg on-ART)."
  ),
  comp4_toxo_pos_vs_2026neg_CC = list(
    keys0 = NEG2026_CC, keys1 = POS_CC,
    label0 = "2026neg_offART_CC", label1 = "pos_offART_CC",
    note   = "Toxo contrast (CC), no confound."
  )
)

# ---- core diff runner ----------------------------------------------------------
run_diff <- function(keys, treatment, covariates, mpg = 3L) {
  files <- as.list(file.path(MKDIR, paste0(keys, "_methylkit.txt")))
  ids   <- as.list(keys)
  obj <- methRead(location = files, sample.id = ids, assembly = "hg38",
                  treatment = treatment, context = "CpG", mincov = 10)
  f  <- filterByCoverage(obj, lo.count = 10, hi.perc = 99.9)
  n  <- normalizeCoverage(f)
  ut <- unite(tileMethylCounts(n), min.per.group = as.integer(mpg))
  calculateDiffMeth(ut, overdispersion = "MN", covariates = covariates, mc.cores = 8)
}

# build covariate data.frame for a set of keys + a model spec
cov_for <- function(keys, cols) {
  if (is.null(cols)) return(NULL)
  p <- prop[match(keys, prop$key)]
  as.data.frame(p[, ..cols])
}

run_comparison <- function(name, def) {
  cat("\n==================================================\n")
  cat("COMPARISON:", name, "\n")
  cat("  0 =", def$label0, "(n=", length(def$keys0), ")\n")
  cat("  1 =", def$label1, "(n=", length(def$keys1), ")\n")
  cat("  NOTE:", def$note, "\n")
  cat("==================================================\n")
  keys <- c(def$keys0, def$keys1)
  treatment <- c(rep(0, length(def$keys0)), rep(1, length(def$keys1)))

  # covariate-vs-status correlations (diagnostic of entanglement)
  cat("\nCovariate-vs-status correlations:\n")
  for (cl in c("Oligo_NNLS","Micro_NNLS","Neuron_total_NNLS")) {
    v <- prop[match(keys, prop$key)][[cl]]
    cat(sprintf("  %-18s r = % .3f\n", cl, cor(v, treatment)))
  }

  # four models, full tables, no cutoff
  for (mname in names(MODELS)) {
    cov <- cov_for(keys, MODELS[[mname]])
    dm  <- run_diff(keys, treatment, cov)
    d   <- getData(dm)
    fwrite(d, file.path(OUTDIR, paste0(name, "__", mname, "__all_tiles.csv")))
    cat(sprintf("  %-13s -> %d tiles  (q<0.05: %d, q<0.10: %d)\n",
                mname, nrow(d), sum(d$qvalue < 0.05), sum(d$qvalue < 0.10)))
  }

  # LOO on all four models
  for (mname in names(MODELS)) {
    cat("\n  LOO on", mname, ":\n")
    loo_rows <- list()
    for (i in seq_along(keys)) {
      keep <- setdiff(seq_along(keys), i)
      kk <- keys[keep]; tt <- treatment[keep]
      n0 <- sum(tt == 0); n1 <- sum(tt == 1)
      mpg <- min(3L, n0, n1)
      cov <- cov_for(kk, MODELS[[mname]])
      dm <- tryCatch(run_diff(kk, tt, cov, mpg = mpg),
                     error = function(e) { cat("    drop", keys[i], "ERR:", conditionMessage(e), "\n"); NULL })
      if (!is.null(dm)) {
        d <- getData(dm)
        loo_rows[[keys[i]]] <- data.table(dropped = keys[i], mpg = mpg,
                                          n_tiles = nrow(d),
                                          n_q05 = sum(d$qvalue < 0.05),
                                          n_q10 = sum(d$qvalue < 0.10))
        fwrite(d, file.path(QCDIR, paste0(name, "__LOO_", mname, "__drop_", keys[i], "__all_tiles.csv")))
      }
    }
    loo_tab <- rbindlist(loo_rows)
    print(loo_tab)
    fwrite(loo_tab, file.path(QCDIR, paste0(name, "__LOO_", mname, "__summary.csv")))
  }
}

for (nm in names(comparisons)) run_comparison(nm, comparisons[[nm]])
cat("\nAll comparisons complete. Full tables in", OUTDIR, "; LOO in", QCDIR, "\n")
sink(type = "message"); sink(); close(log_con)
