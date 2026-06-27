suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

# ── Covariates from scMD proportions ─────────────────────────────────────────
prop <- fread("scMD_proportions.tsv")
setnames(prop, 1, "sample")
sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
prop <- prop[match(sample_order, prop$sample)]
stopifnot(all(prop$sample == sample_order))

exc        <- prop$Exc
oligo_lin  <- prop$Oligo + prop$OPC
micro      <- prop$Micro

treatment <- c(0,0,0,0,0,0,1,1,1)   # 0=toxo+, 1=toxo-

# Report each covariate's correlation with treatment (collinearity watch)
cat("Covariate correlations with treatment:\n")
cat("  Exc:          ", round(cor(exc, treatment),3), "\n")
cat("  Oligo-lineage:", round(cor(oligo_lin, treatment),3), "\n")
cat("  Micro:        ", round(cor(micro, treatment),3), "\n\n")

# ── methylKit base pipeline (once) ───────────────────────────────────────────
sample_files <- as.list(paste0(sample_order, "_methylkit.txt"))
sample_ids   <- as.list(sample_order)
myobj <- methRead(location=sample_files, sample.id=sample_ids,
                  assembly="hg38", treatment=treatment, context="CpG", mincov=10)
filtered     <- filterByCoverage(myobj, lo.count=10, hi.perc=99.9)
normalized   <- normalizeCoverage(filtered)
united_tiles <- tileMethylCounts(normalized) |> unite(min.per.group=3L)

# ── Define the models ────────────────────────────────────────────────────────
models <- list(
    unadjusted        = NULL,
    exc_only          = data.frame(exc=exc),
    oligolin_only     = data.frame(oligo_lineage=oligo_lin),
    micro_only        = data.frame(micro=micro),
    exc_micro         = data.frame(exc=exc, micro=micro),
    exc_oligolin      = data.frame(exc=exc, oligo_lineage=oligo_lin)
)

# Reference set = unadjusted significant tiles (for survived/lost/new)
run_model <- function(cov) {
    if (is.null(cov)) calculateDiffMeth(united_tiles, overdispersion="MN", mc.cores=8)
    else calculateDiffMeth(united_tiles, overdispersion="MN", covariates=cov, mc.cores=8)
}

key <- function(x) paste(getData(x)$chr, getData(x)$start, getData(x)$end)

cat("Running", length(models), "models...\n\n")
dm <- lapply(models, run_model)
sig <- lapply(dm, function(d) getMethylDiff(d, difference=10, qvalue=0.05))

ref_key <- key(sig$unadjusted)

# ── Summary table ────────────────────────────────────────────────────────────
cat(sprintf("%-16s %6s %9s %6s %5s\n", "model", "n_sig", "survived", "lost", "new"))
cat(strrep("-", 50), "\n")
for (nm in names(models)) {
    k <- key(sig[[nm]])
    if (nm == "unadjusted") {
        cat(sprintf("%-16s %6d %9s %6s %5s\n", nm, length(k), "-", "-", "-"))
    } else {
        cat(sprintf("%-16s %6d %9d %6d %5d\n", nm, length(k),
            sum(ref_key %in% k), sum(!ref_key %in% k), sum(!k %in% ref_key)))
    }
}

# ── Write outputs ────────────────────────────────────────────────────────────
for (nm in names(models)) {
    fwrite(getData(dm[[nm]]),  paste0("all_tiles_", nm, ".csv"))      # all ~86k
    fwrite(getData(sig[[nm]]), paste0("sig_", nm, ".csv"))           # q<0.05 subset
}
cat("\nWrote all_tiles_*.csv (full) and sig_*.csv (q<0.05) for each model.\n")

# Cross-model survival of the unadjusted 293: which models keep each tile
ref_data <- getData(sig$unadjusted)
surv <- data.table(chr=ref_data$chr, start=ref_data$start, end=ref_data$end,
                   meth.diff=ref_data$meth.diff, qvalue=ref_data$qvalue)
for (nm in setdiff(names(models),"unadjusted")) {
    surv[[nm]] <- key(sig$unadjusted) %in% key(sig[[nm]])
}
surv[, n_models_surviving := rowSums(.SD), .SDcols=setdiff(names(models),"unadjusted")]
fwrite(surv, "unadj_tile_survival_across_models.csv")
cat("Wrote unadj_tile_survival_across_models.csv (each of 293 tiles x which models retain it)\n")
