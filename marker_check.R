suppressPackageStartupMessages(library(scMD))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

sig <- Lee_sig_all_WGBS

sig_dt <- data.table(coord_raw = rownames(sig))
sig_dt[, chr := paste0("chr", tstrsplit(coord_raw, ":")[[1]])]
sig_dt[, pos := as.integer(tstrsplit(coord_raw, ":")[[2]])]
sig_dt[, coord := paste(chr, pos, sep=":")]
all_ct <- colnames(sig)
for (ct in all_ct) sig_dt[[ct]] <- sig[, ct]

# ── cell-type specificity: target far from mean of the others ────────────────
add_delta <- function(dt, target, others) {
    om <- rowMeans(as.matrix(dt[, ..others]))
    dt[[target]] - om
}
sig_dt[, astro_delta := add_delta(sig_dt, "Astro", setdiff(all_ct,"Astro"))]
sig_dt[, micro_delta := add_delta(sig_dt, "Micro", setdiff(all_ct,"Micro"))]

# ── load all 9 samples, keyed by coord ───────────────────────────────────────
sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
bed_list <- lapply(sample_order, function(s) {
    b <- fread(file.path(s, paste0(s,"_GRCh38_CpG.bed")),
               col.names=c("chr","start","end","beta","coverage"))
    b[, coord := paste(chr, start, sep=":")]; setkey(b, coord); b
})
names(bed_list) <- sample_order

# ── for a set of candidate coords, return per-sample beta + coverage ──────────
eval_markers <- function(cand, label) {
    cov  <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), coverage])
    beta <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), beta])
    rownames(cov) <- cand$coord; rownames(beta) <- cand$coord
    universal <- which(rowSums(!is.na(cov)) == length(sample_order))
    usable <- integer(0)
    for (i in universal) if (sd(beta[i,])>0.1 && all(cov[i,]>=10)) usable <- c(usable, i)
    cat(sprintf("\n%s: %d specific(|delta|>0.3), %d universal, %d usable(depth>=10 & range)\n",
                label, nrow(cand), length(universal), length(usable)))
    if (length(usable)>0) {
        u <- cand$coord[usable]
        cat("Usable marker coords:\n"); print(u)
        cat("Per-sample depth:\n");    print(cov[u, , drop=FALSE])
        cat("Per-sample beta:\n");     print(round(beta[u, , drop=FALSE], 3))
    }
}

astro_cand <- sig_dt[abs(astro_delta) > 0.6]
micro_cand <- sig_dt[abs(micro_delta) > 0.6]

eval_markers(astro_cand, "Astrocyte")
eval_markers(micro_cand, "Microglia")
