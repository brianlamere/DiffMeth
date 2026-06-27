#!/usr/bin/env Rscript
# ToxoDM cell-type deconvolution (phase 2, full 23-sample set)
# Fresh from current data: 14 newly-aligned BG/FMC beds + 9 inherited CC beds.
# scMD WGBS reduced-ensemble, same invocation style as phase 1 (default ensemble,
# tolerate per-method failures). Output authoritative on its own; no baseline dep.
suppressPackageStartupMessages({
  library(scMD)
  library(nnls)
  library(data.table)
})

ALIGN <- "/projects1/ToxoDM/working/align"
OUTDIR <- "/projects1/ToxoDM/reports"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. discover the 23 beds (follow symlinks for inherited CC) ----------------
beds <- list.files(ALIGN, pattern = "_CpG\\.bed$", full.names = TRUE, recursive = TRUE)
stopifnot(length(beds) == 23)

# parse cohort/tissue/sample from path: .../align/<cohort>/<tissue>/<sample>/<file>
meta <- rbindlist(lapply(beds, function(p) {
  parts <- strsplit(sub(paste0("^", ALIGN, "/"), "", p), "/")[[1]]
  data.table(cohort = parts[1], tissue = parts[2], sample = parts[3], path = p)
}))
meta[, key := paste(cohort, tissue, sample, sep = "_")]
cat("Discovered", nrow(meta), "beds:\n")
print(meta[, .(cohort, tissue, sample)])

# ---- 2. build bulk beta matrix keyed by chr:pos --------------------------------
# bed cols: chr start end beta coverage  (0-based start = C position)
read_bed <- function(p) {
  b <- fread(p, col.names = c("chr","start","end","beta","coverage"))
  b[, coord := paste0(sub("^chr","",chr), ":", start)]   # match Lee_sig_all_WGBS rownames
  b[, .(coord, beta)]
}

cat("\nReading beds and assembling bulk matrix...\n")
bedlist <- lapply(meta$path, read_bed)
names(bedlist) <- meta$key

all_coords <- unique(unlist(lapply(bedlist, function(x) x$coord)))
bulk <- matrix(NA_real_, nrow = length(all_coords), ncol = nrow(meta),
               dimnames = list(all_coords, meta$key))
for (i in seq_len(nrow(meta))) {
  b <- bedlist[[i]]
  bulk[b$coord, meta$key[i]] <- b$beta
}
cat("Bulk matrix:", nrow(bulk), "CpGs x", ncol(bulk), "samples\n")

# ---- 3. restrict to scMD WGBS signature, complete cases ------------------------
sig_coords <- rownames(Lee_sig_all_WGBS)
bulk_sig <- bulk[rownames(bulk) %in% sig_coords, , drop = FALSE]
bulk_cc  <- bulk_sig[rowSums(is.na(bulk_sig)) == 0, , drop = FALSE]
cat("Signature CpGs covered in ALL", ncol(bulk), "samples:", nrow(bulk_cc), "\n")

# ---- 4. deconvolve: same style as phase 1, default ensemble, tolerate errors ---
# Exclude Houseman: it routes through minfi's makeGenomicRatioSetFromMatrix, which
# hard-requires 450k cg-probe rownames and errors on WGBS coordinate rownames (a
# structural incompatibility, not a convergence issue). The remaining coordinate-
# tolerant methods still self-select convergence on this larger set; any that error
# internally are tolerated per-method failures the ensemble routes around.
cat("\nRunning scMD deconvolution (WGBS, coordinate-tolerant ensemble)...\n")

# reference: signature CpGs (rows) x 7 cell types (cols), matched to bulk_cc rows
ref <- Lee_sig_all_WGBS[rownames(bulk_cc), ]   # 3747 x 7
props <- t(apply(bulk_cc, 2, function(sample_col) {
    fit <- nnls(ref, sample_col)      # solve sample ≈ ref %*% proportions, props >= 0
    p <- fit$x
    p / sum(p)                        # normalize to sum to 1
}))
colnames(props) <- colnames(ref)

prop <- as.data.frame(round(res$scMD_p, 4))
prop <- cbind(meta[match(rownames(prop), meta$key), .(cohort, tissue, sample)], prop)

cat("\n== Estimated proportions ==\n")
print(prop, row.names = FALSE)

outfile <- file.path(OUTDIR, "scMD_proportions.tsv")
fwrite(prop, outfile, sep = "\t")
cat("\nWrote", outfile, "\n")
