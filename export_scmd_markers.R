suppressPackageStartupMessages(library(scMD))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

sig <- Lee_sig_all_WGBS
sig_dt <- data.table(coord_raw = rownames(sig))
sig_dt[, chr := paste0("chr", tstrsplit(coord_raw, ":")[[1]])]
sig_dt[, pos := as.integer(tstrsplit(coord_raw, ":")[[2]])]
sig_dt[, coord := paste(chr, pos, sep=":")]
all_ct <- colnames(sig)
for (ct in all_ct) set(sig_dt, j=ct, value=sig[, ct])

# delta (specificity) per target cell type
for (target in c("Astro","Micro","Oligo")) {
    others <- setdiff(all_ct, target)
    om <- rowMeans(as.matrix(sig_dt[, ..others]))
    set(sig_dt, j=paste0(target,"_delta"), value=sig_dt[[target]] - om)
}

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
cohort       <- c(rep("2024",6), rep("2026",3))
bed_list <- lapply(sample_order, function(s) {
    b <- fread(file.path(s, paste0(s,"_GRCh38_CpG.bed")),
               col.names=c("chr","start","end","beta","coverage"))
    b[, coord := paste(chr, start, sep=":")]; setkey(b, coord); b
})
names(bed_list) <- sample_order

DELTA <- 0.5

# For one cell type: markers specific at DELTA, covered in all 9 samples.
# Returns long-format rows: coord, celltype, ref betas, per-sample beta+depth.
export_celltype <- function(target) {
    cand <- sig_dt[abs(get(paste0(target,"_delta"))) > DELTA]
    cov  <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), coverage])
    beta <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), beta])
    keep <- which(rowSums(!is.na(cov)) == length(sample_order))   # covered in all 9
    if (length(keep)==0) return(NULL)
    out <- data.table(
        celltype       = target,
        coord          = cand$coord[keep],
        chr            = cand$chr[keep],
        pos            = cand$pos[keep],
        ref_target     = round(cand[[target]][keep], 4),
        target_delta   = round(cand[[paste0(target,"_delta")]][keep], 4)
    )
    # also flag whether this marker is ALSO specific for the other two glial types
    for (other in setdiff(c("Astro","Micro","Oligo"), target))
        out[[paste0("also_",other)]] <- abs(cand[[paste0(other,"_delta")]][keep]) > DELTA
    # per-sample beta and depth columns
    bk <- beta[keep, , drop=FALSE]; ck <- cov[keep, , drop=FALSE]
    for (i in seq_along(sample_order)) {
        out[[paste0(sample_order[i], "_beta")]]  <- round(bk[, i], 4)
        out[[paste0(sample_order[i], "_depth")]] <- ck[, i]
    }
    out
}

astro <- export_celltype("Astro")
micro <- export_celltype("Micro")
oligo <- export_celltype("Oligo")

combined <- rbindlist(list(astro, micro, oligo), fill=TRUE)
cat("Markers exported (covered in all 9, |delta|>", DELTA, "):\n", sep="")
cat("  Astro:", nrow(astro), "  Micro:", nrow(micro), "  Oligo:", nrow(oligo), "\n")

# ── write CSVs ────────────────────────────────────────────────────────────────
# 1. The combined long table: one row per marker per cell type, ref + per-sample beta & depth
fwrite(combined, "scMD_markers_full.csv")

# 2. A tidy beta-only matrix (markers x samples) for quick stats in a spreadsheet
beta_wide <- combined[, c("celltype","coord","chr","pos",
    paste0(sample_order,"_beta")), with=FALSE]
fwrite(beta_wide, "scMD_markers_beta.csv")

# 3. A tidy depth-only matrix
depth_wide <- combined[, c("celltype","coord","chr","pos",
    paste0(sample_order,"_depth")), with=FALSE]
fwrite(depth_wide, "scMD_markers_depth.csv")

cat("\nWrote:\n  scMD_markers_full.csv  (ref + per-sample beta & depth, with cross-celltype flags)\n")
cat("  scMD_markers_beta.csv  (markers x sample beta matrix)\n")
cat("  scMD_markers_depth.csv (markers x sample depth matrix)\n")
