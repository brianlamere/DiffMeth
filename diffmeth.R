suppressPackageStartupMessages(library(methylKit))

setwd("/projects2/DiffMeth")

# ── 1. Glial fractions from HiBED deconvolution (deconv_buckets.R) ───────────
# Bucket1_2024: CA192,CA346,CB239,CC249,CE167,CE234 — Glial 0.8491
# Bucket2_2026: LG30,LG31,LG52                      — Glial 0.5698
glial_fractions <- data.frame(
    glial = c(0.8491, 0.8491, 0.8491, 0.8491, 0.8491, 0.8491,
              0.5698, 0.5698, 0.5698)
)

# ── 2. Load methylKit txt files ───────────────────────────────────────────────
sample_files <- list(
    "CA192_methylkit.txt", "CA346_methylkit.txt", "CB239_methylkit.txt",
    "CC249_methylkit.txt", "CE167_methylkit.txt", "CE234_methylkit.txt",
    "LG30_methylkit.txt",  "LG31_methylkit.txt",  "LG52_methylkit.txt"
)

sample_ids <- list(
    "CA192","CA346","CB239","CC249","CE167","CE234",
    "LG30","LG31","LG52"
)

treatment <- c(0,0,0,0,0,0,1,1,1)

myobj <- methRead(
    location   = sample_files,
    sample.id  = sample_ids,
    assembly   = "hg38",
    treatment  = treatment,
    context    = "CpG",
    mincov     = 10
)

# ── 3. Filter, normalize, unite ───────────────────────────────────────────────
filtered   <- filterByCoverage(myobj, lo.count=10, hi.perc=99.9)

normalized <- normalizeCoverage(filtered)
united     <- unite(normalized, min.per.group=3L)
united_tiles <- tileMethylCounts(normalized) |>
    unite(min.per.group=3L)

# ── 4. adding QC:  answering Sarah's question ───────────────────────────────────────────

# Q1: filterByCoverage parameters used
cat("\n── Coverage filter parameters ──\n")
cat("lo.count = 10 (minimum coverage)\n")
cat("hi.perc = 99.9 (top 0.1% excluded as potential PCR artifacts)\n")

# Q2: Coverage distribution before and after filtering
pdf("methylkit_QC_coverage.pdf", width=12, height=8)

# Before filtering - per sample coverage summary
cat("\n── Coverage summary BEFORE filtering ──\n")
for (i in seq_along(myobj)) {
    d <- getData(myobj[[i]])
    cat(sample_ids[[i]], "- median:", median(d$coverage),
        "mean:", round(mean(d$coverage), 1),
        "max:", max(d$coverage), "\n")
}

# Plot coverage distributions before filtering
par(mfrow=c(3,3), mar=c(4,4,2,1))
for (i in seq_along(myobj)) {
    d <- getData(myobj[[i]])
    hist(log10(d$coverage + 1), breaks=50,
         main=paste(sample_ids[[i]], "pre-filter"),
         xlab="log10(coverage + 1)", col="steelblue", border=NA)
}
dev.off()

# After filtering
pdf("methylkit_QC_coverage_postfilter.pdf", width=12, height=8)
cat("\n── Coverage summary AFTER filtering ──\n")
for (i in seq_along(filtered)) {
    d <- getData(filtered[[i]])
    cat(sample_ids[[i]], "- median:", median(d$coverage),
        "mean:", round(mean(d$coverage), 1),
        "max:", max(d$coverage),
        "sites retained:", nrow(d), "\n")
}

par(mfrow=c(3,3), mar=c(4,4,2,1))
for (i in seq_along(filtered)) {
    d <- getData(filtered[[i]])
    hist(log10(d$coverage + 1), breaks=50,
         main=paste(sample_ids[[i]], "post-filter"),
         xlab="log10(coverage + 1)", col="coral", border=NA)
}
dev.off()

# Q3: Per-tile coverage by sample from united_tiles
cat("\n── Per-tile coverage summary (united_tiles) ──\n")
tile_data <- getData(united_tiles)
cov_cols  <- grep("coverage", names(tile_data), value=TRUE)
cat("Coverage columns:", paste(cov_cols, collapse=", "), "\n")
print(summary(tile_data[, cov_cols]))

# Write per-tile coverage table
cov_table <- tile_data[, c("chr","start","end", cov_cols)]
write.table(cov_table, "methylkit_tile_coverage_by_sample.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Per-tile coverage written to methylkit_tile_coverage_by_sample.tsv\n")

# Normalization check - global methylation before/after
pdf("methylkit_QC_normalization.pdf", width=10, height=5)
par(mfrow=c(1,2))

# Before normalization
pre_meth <- sapply(filtered, function(x) {
    d <- getData(x)
    sum(d$numCs) / sum(d$coverage)
})
names(pre_meth) <- unlist(sample_ids)
barplot(pre_meth * 100, las=2,
        col=c(rep("steelblue",6), rep("coral",3)),
        ylab="Global CpG methylation %",
        main="Before normalization",
        ylim=c(0, max(pre_meth * 100) * 1.2))

# After normalization
post_meth <- sapply(normalized, function(x) {
    d <- getData(x)
    sum(d$numCs) / sum(d$coverage)
})
names(post_meth) <- unlist(sample_ids)
barplot(post_meth * 100, las=2,
        col=c(rep("steelblue",6), rep("coral",3)),
        ylab="Global CpG methylation %",
        main="After normalization",
        ylim=c(0, max(post_meth * 100) * 1.2))
dev.off()

cat("\nQC files written:\n")
cat("  methylkit_QC_coverage.pdf\n")
cat("  methylkit_QC_coverage_postfilter.pdf\n")
cat("  methylkit_QC_normalization.pdf\n")
cat("  methylkit_tile_coverage_by_sample.tsv\n")

# ── 5. Unadjusted differential methylation ───────────────────────────────────
message("Running unadjusted DM...")
dm_tiles_unadj <- calculateDiffMeth(
    united_tiles,
    overdispersion = "MN",
    mc.cores       = 8
)


#sig_unadj_q010       <- getMethylDiff(dm_tiles_unadj, difference=10, qvalue=0.10)
#write.table(getData(sig_unadj_q010), "sig_tiles_unadjusted_q010.tsv",
#            sep="\t", quote=FALSE, row.names=FALSE)

sig_unadj       <- getMethylDiff(dm_tiles_unadj, difference=10, qvalue=0.05)
sig_unadj_hyper <- getMethylDiff(dm_tiles_unadj, difference=10, qvalue=0.05, type="hyper")
sig_unadj_hypo  <- getMethylDiff(dm_tiles_unadj, difference=10, qvalue=0.05, type="hypo")

cat("\n── Unadjusted results ──\n")
cat("Total significant tiles:", nrow(sig_unadj), "\n")
cat("Hyper (2026 > 2024):", nrow(sig_unadj_hyper), "\n")
cat("Hypo  (2026 < 2024):", nrow(sig_unadj_hypo), "\n")


# ── 6. Glial-adjusted differential methylation ───────────────────────────────
message("Running glial-adjusted DM...")
dm_tiles_adj <- calculateDiffMeth(
    united_tiles,
    overdispersion = "MN",
    covariates     = glial_fractions,
    mc.cores       = 8
)

sig_adj       <- getMethylDiff(dm_tiles_adj, difference=10, qvalue=0.05)
sig_adj_hyper <- getMethylDiff(dm_tiles_adj, difference=10, qvalue=0.05, type="hyper")
sig_adj_hypo  <- getMethylDiff(dm_tiles_adj, difference=10, qvalue=0.05, type="hypo")

cat("\n── Glial-adjusted results ──\n")
cat("Total significant tiles:", nrow(sig_adj), "\n")
cat("Hyper (2026 > 2024):", nrow(sig_adj_hyper), "\n")
cat("Hypo  (2026 < 2024):", nrow(sig_adj_hypo), "\n")

# ── 7. Comparison: which tiles survive adjustment ─────────────────────────────
unadj_coords <- paste(getData(sig_unadj)$chr, getData(sig_unadj)$start)
adj_coords   <- paste(getData(sig_adj)$chr,   getData(sig_adj)$start)

cat("\n── Comparison ──\n")
cat("Tiles surviving adjustment:", sum(unadj_coords %in% adj_coords), "\n")
cat("Tiles lost after adjustment:", sum(!unadj_coords %in% adj_coords), "\n")
cat("New tiles appearing after adjustment:", sum(!adj_coords %in% unadj_coords), "\n")

# ── 8. Volcano-style comparison plot ─────────────────────────────────────────
unadj_data <- getData(dm_tiles_unadj)
adj_data   <- getData(dm_tiles_adj)

merged <- merge(
    unadj_data[, c("chr","start","end","meth.diff","qvalue")],
    adj_data[,   c("chr","start","end","meth.diff","qvalue")],
    by     = c("chr","start","end"),
    suffixes = c("_unadj","_adj")
)

pdf("diffmeth_comparison.pdf", width=8, height=8)

# Unadjusted vs adjusted meth.diff scatter
plot(
    merged$meth.diff_unadj,
    merged$meth.diff_adj,
    pch  = 20,
    cex  = 0.2,
    col  = "grey60",
    xlab = "Methylation difference unadjusted (2026 - 2024)",
    ylab = "Methylation difference adjusted for glial fraction",
    main = "Effect of glial fraction adjustment on DM estimates"
)
abline(0, 1, lty=2, col="red")

# Highlight tiles significant in both
both_sig <- merged$qvalue_unadj < 0.05 & merged$meth.diff_unadj > 10 |
            merged$qvalue_unadj < 0.05 & merged$meth.diff_unadj < -10
points(
    merged$meth.diff_unadj[both_sig],
    merged$meth.diff_adj[both_sig],
    pch = 20, cex = 0.5, col = "red"
)

dev.off()

# ── 9. Save results tables ────────────────────────────────────────────────────
write.table(getData(sig_unadj), "sig_tiles_unadjusted.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(getData(sig_adj),   "sig_tiles_adjusted.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(merged, "dm_tiles_unadj_vs_adj.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

message("Done. Output files written to /projects2/DiffMeth/")
