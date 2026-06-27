#!/usr/bin/env Rscript
# ToxoDM phase-2 differential methylation battery — 3 comparisons.
# Per comparison: unadjusted + primary(neuron+oligo-lineage) + sensitivity variants
# + LOO. Covariates from scMD_proportions.tsv. All q-values reported, NO cutoff.
# Self-contained: reads methylKit txt + proportions from disk, no session state.
suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
})

MKDIR  <- "/projects1/ToxoDM/working/methylkit"
PROPF  <- "/projects1/ToxoDM/reports/scMD_proportions.tsv"
OUTDIR <- "/projects1/ToxoDM/reports"
QCDIR  <- file.path(OUTDIR, "QC")
dir.create(QCDIR, showWarnings = FALSE, recursive = TRUE)

prop <- fread(PROPF)
prop[, key := paste(cohort, tissue, sample, sep = "_")]

# Covariate builders from proportions
covsets <- function(keys) {
  p <- prop[match(keys, prop$key)]
  list(
    none          = NULL,
    primary       = data.frame(neuron = p$Inh + p$Exc, oligolin = p$Oligo + p$OPC),
    exc_only      = data.frame(exc = p$Exc),
    oligolin_only = data.frame(oligolin = p$Oligo + p$OPC),
    micro_only    = data.frame(micro = p$Micro),
    exc_micro     = data.frame(exc = p$Exc, micro = p$Micro),
    exc_oligolin  = data.frame(exc = p$Exc, oligolin = p$Oligo + p$OPC)
  )
}

# ---- comparison definitions ----------------------------------------------------
# treatment polarity: reference group = 0, treatment group = 1
# positives = 6 distinct donors (CE234 read for the email's duplicate-CB239 typo)
NEG_FMC <- c("2024_FMC_Atpsy_15","2024_FMC_Atpsy_20","2024_FMC_Atpsy_21","2024_FMC_LG29")
NEG_BG  <- c("2024_BG_Atpsy_15","2024_BG_Atpsy_20","2024_BG_Atpsy_21","2024_BG_LG29")
POS_CC  <- c("2024_CC_CA192","2024_CC_CA346","2024_CC_CB239","2024_CC_CC249","2024_CC_CE167","2024_CC_CE234")
POS_BG  <- c("2024_BG_CA192","2024_BG_CA346","2024_BG_CB239","2024_BG_CC249","2024_BG_CE167","2024_BG_CE234")
NEG_CC2026 <- c("2026_CC_LG30","2026_CC_LG31","2026_CC_LG52")

comparisons <- list(
  comp1_neg2024_vs_neg2026 = list(
    keys0 = NEG_FMC,    keys1 = NEG_CC2026,
    label0 = "2024neg_FMC", label1 = "2026neg_CC"
  ),
  comp2_pos_vs_neg_CC = list(
    keys0 = NEG_FMC,    keys1 = POS_CC,
    label0 = "neg_FMC", label1 = "pos_CC"
  ),
  comp3_pos_vs_neg_BG = list(
    keys0 = NEG_BG,     keys1 = POS_BG,
    label0 = "neg_BG",  label1 = "pos_BG"
  )
)

# ---- core runner ---------------------------------------------------------------
run_diff <- function(keys, treatment, covariates, mpg = 3L) {
  files <- as.list(file.path(MKDIR, paste0(keys, "_methylkit.txt")))
  ids   <- as.list(keys)
  obj <- methRead(location = files, sample.id = ids, assembly = "hg38",
                  treatment = treatment, context = "CpG", mincov = 10)
  f <- filterByCoverage(obj, lo.count = 10, hi.perc = 99.9)
  n <- normalizeCoverage(f)
  ut <- unite(tileMethylCounts(n), min.per.group = as.integer(mpg))
  calculateDiffMeth(ut, overdispersion = "MN",
                    covariates = covariates, mc.cores = 8)
}

run_comparison <- function(name, def) {
  cat("\n==================================================\n")
  cat("COMPARISON:", name, "\n")
  cat("  0 =", def$label0, "(n=", length(def$keys0), ")\n")
  cat("  1 =", def$label1, "(n=", length(def$keys1), ")\n")
  cat("==================================================\n")
  keys <- c(def$keys0, def$keys1)
  treatment <- c(rep(0, length(def$keys0)), rep(1, length(def$keys1)))
  cov <- covsets(keys)

  # collinearity diagnostic: each covariate column vs treatment
  cat("\nCovariate-vs-treatment correlations:\n")
  for (nm in setdiff(names(cov), "none")) {
    cm <- cov[[nm]]
    for (cc in colnames(cm)) {
      cat(sprintf("  %-14s %-10s r = % .3f\n", nm, cc, cor(cm[[cc]], treatment)))
    }
  }

  # run all models, dump full results (no cutoff)
  for (nm in names(cov)) {
    dm <- run_diff(keys, treatment, cov[[nm]])
    fwrite(getData(dm), file.path(OUTDIR, paste0(name, "__", nm, "__all_tiles.csv")))
    cat(sprintf("  model %-14s -> %d tiles (full, no cutoff)\n", nm, nrow(getData(dm))))
  }

  # LOO on unadjusted + primary
  for (model in c("none","primary")) {
    cat("\n  LOO on", model, ":\n")
    loo_rows <- list()
    for (i in seq_along(keys)) {
      keep <- setdiff(seq_along(keys), i)
      kk <- keys[keep]; tt <- treatment[keep]
      n0 <- sum(tt == 0); n1 <- sum(tt == 1)
      mpg <- min(3L, n0, n1)
      cset <- if (model == "none") NULL else covsets(kk)[["primary"]]
      dm <- tryCatch(run_diff(kk, tt, cset, mpg = mpg),
                     error = function(e) { cat("    drop", keys[i], "ERR:", conditionMessage(e), "\n"); NULL })
      if (!is.null(dm)) {
        d <- getData(dm)
        loo_rows[[keys[i]]] <- data.table(dropped = keys[i], mpg = mpg,
                                          n_tiles = nrow(d),
                                          n_q05 = sum(d$qvalue < 0.05),
                                          n_q10 = sum(d$qvalue < 0.10))
        fwrite(d, file.path(QCDIR, paste0(name, "__LOO_", model, "__drop_", keys[i], "__all_tiles.csv")))
      }
    }
    loo_tab <- rbindlist(loo_rows)
    print(loo_tab)
    fwrite(loo_tab, file.path(QCDIR, paste0(name, "__LOO_", model, "__summary.csv")))
  }
}

for (nm in names(comparisons)) run_comparison(nm, comparisons[[nm]])
cat("\nAll comparisons complete. Full tables in", OUTDIR, "; LOO in", QCDIR, "\n")
