suppressPackageStartupMessages(library(scMD))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

sig <- Lee_sig_all_WGBS
sig_dt <- data.table(coord_raw = rownames(sig))
sig_dt[, chr := paste0("chr", tstrsplit(coord_raw, ":")[[1]])]
sig_dt[, pos := as.integer(tstrsplit(coord_raw, ":")[[2]])]
sig_dt[, coord := paste(chr, pos, sep=":")]
all_ct <- colnames(sig)   # Astro Micro Endo Oligo OPC Inh Exc
for (ct in all_ct) set(sig_dt, j=ct, value=sig[, ct])

# specificity delta for every cell type
for (target in all_ct) {
    others <- setdiff(all_ct, target)
    om <- rowMeans(as.matrix(sig_dt[, ..others]))
    set(sig_dt, j=paste0(target,"_delta"), value=sig_dt[[target]] - om)
}

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
bed_list <- lapply(sample_order, function(s) {
    b <- fread(file.path(s, paste0(s,"_GRCh38_CpG.bed")),
               col.names=c("chr","start","end","beta","coverage"))
    b[, coord := paste(chr, start, sep=":")]; setkey(b, coord); b
})
names(bed_list) <- sample_order

DELTA <- 0.5

export_celltype <- function(target) {
    cand <- sig_dt[abs(get(paste0(target,"_delta"))) > DELTA]
    cov  <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), coverage])
    beta <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), beta])
    keep <- which(rowSums(!is.na(cov)) == length(sample_order))
    if (length(keep)==0) return(NULL)
    out <- data.table(
        celltype     = target,
        coord        = cand$coord[keep],
        chr          = cand$chr[keep],
        pos          = cand$pos[keep],
        target_delta = round(cand[[paste0(target,"_delta")]][keep], 4)
    )
    # reference betas for ALL 7 cell types
    for (ct in all_ct)
        out[[paste0("ref_", ct)]] <- round(cand[[ct]][keep], 4)
    # flags: also specific for each OTHER cell type at this delta
    for (ct in setdiff(all_ct, target))
        out[[paste0("also_", ct)]] <- abs(cand[[paste0(ct,"_delta")]][keep]) > DELTA
    # per-sample beta and depth
    bk <- beta[keep, , drop=FALSE]; ck <- cov[keep, , drop=FALSE]
    for (i in seq_along(sample_order)) {
        out[[paste0(sample_order[i], "_beta")]]  <- round(bk[, i], 4)
        out[[paste0(sample_order[i], "_depth")]] <- ck[, i]
    }
    out
}

per_ct <- lapply(all_ct, export_celltype)
names(per_ct) <- all_ct
combined <- rbindlist(per_ct, fill=TRUE)

cat("Markers per cell type (|delta|>", DELTA, ", covered in all 9):\n", sep="")
for (ct in all_ct) cat(sprintf("  %-6s %d\n", ct, nrow(per_ct[[ct]])))
cat("  TOTAL rows (markers may repeat across types):", nrow(combined), "\n")
cat("  unique marker coords:", length(unique(combined$coord)), "\n")

fwrite(combined, "scMD_markers_full_7ct.csv")

beta_wide  <- combined[, c("celltype","coord","chr","pos", paste0(sample_order,"_beta")),  with=FALSE]
depth_wide <- combined[, c("celltype","coord","chr","pos", paste0(sample_order,"_depth")), with=FALSE]
fwrite(beta_wide,  "scMD_markers_beta_7ct.csv")
fwrite(depth_wide, "scMD_markers_depth_7ct.csv")

cat("\nWrote:\n")
cat("  scMD_markers_full_7ct.csv  (target_delta + all-7 ref betas + cross-type flags + per-sample beta & depth)\n")
cat("  scMD_markers_beta_7ct.csv  (markers x sample beta)\n")
cat("  scMD_markers_depth_7ct.csv (markers x sample depth)\n")
