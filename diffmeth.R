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

# ── 4. Tile-level objects ─────────────────────────────────────────────────────
united_tiles <- tileMethylCounts(normalized) |>
    unite(min.per.group=3L)

# ── 5. Unadjusted differential methylation ───────────────────────────────────
message("Running unadjusted DM...")
dm_tiles_unadj <- calculateDiffMeth(
    united_tiles,
    overdispersion = "MN",
    mc.cores       = 8
)

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
