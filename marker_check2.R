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

# delta for each target cell type
for (target in c("Astro","Micro","Oligo")) {
    others <- setdiff(all_ct, target)
    om <- rowMeans(as.matrix(sig_dt[, ..others]))
    set(sig_dt, j=paste0(target,"_delta"), value=sig_dt[[target]] - om)
}

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
b2024 <- 1:6; b2026 <- 7:9
bed_list <- lapply(sample_order, function(s) {
    b <- fread(file.path(s, paste0(s,"_GRCh38_CpG.bed")),
               col.names=c("chr","start","end","beta","coverage"))
    b[, coord := paste(chr, start, sep=":")]; setkey(b, coord); b
})
names(bed_list) <- sample_order

DELTA  <- 0.5     # specificity floor
MINDEP <- 10      # depth floor (resolve intermediate values)
MINSD  <- 0.1     # within-sample dynamic range
CLUST  <- 42     # bp window: collapse markers closer than this to one

# candidate markers for a target: specific, NOT specific to either other type
get_candidates <- function(target, exclude) {
    spec <- abs(sig_dt[[paste0(target,"_delta")]]) > DELTA
    # exclude CpGs that are also specific for the other two glial types
    for (ex in exclude) spec <- spec & !(abs(sig_dt[[paste0(ex,"_delta")]]) > DELTA)
    sig_dt[spec, .(coord, chr, pos)]
}

evaluate <- function(cand, label) {
    if (nrow(cand)==0) { cat(label, ": no unique candidates\n"); return(NULL) }
    cov  <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), coverage])
    beta <- sapply(sample_order, function(s) bed_list[[s]][.(cand$coord), beta])
    rownames(cov) <- cand$coord; rownames(beta) <- cand$coord
    keep <- c()
    for (i in seq_len(nrow(cand))) {
        if (any(is.na(cov[i,]))) next
        if (all(cov[i,] >= MINDEP) && sd(beta[i,]) > MINSD) {
            # reject cohort-separated (collinear) markers:
            # require within-cohort variation in BOTH cohorts
            sd24 <- sd(beta[i, b2024]); sd26 <- sd(beta[i, b2026])
            if (sd24 > 0.05 && sd26 > 0.05) keep <- c(keep, i)
        }
    }
    if (length(keep)==0) { cat(label, ": 0 markers pass all filters\n"); return(NULL) }
    cand_keep <- cand[keep]
    # spatial thinning: within each chr, collapse markers within CLUST bp
    cand_keep <- cand_keep[order(chr, pos)]
    cand_keep[, grp := cumsum(c(TRUE, diff(pos) > CLUST | head(chr,-1)!=tail(chr,-1))), by=chr]
    thinned <- cand_keep[, .SD[1], by=.(chr, grp)]
    cat(sprintf("\n%s: %d unique+covered+dynamic, %d after spatial thinning\n",
                label, length(keep), nrow(thinned)))
    print(thinned$coord)
    # build coverage-weighted composite per sample
    co <- thinned$coord
    comp <- sapply(sample_order, function(s) {
        b <- beta[co, s]; w <- cov[co, s]
        sum(b*w)/sum(w)
    })
    cat("Per-sample composite:\n"); print(round(comp, 3))
    cat(sprintf("  2024 mean=%.3f  2026 mean=%.3f  gap=%.3f\n",
        mean(comp[b2024]), mean(comp[b2026]), abs(mean(comp[b2024])-mean(comp[b2026]))))
    list(markers=thinned$coord, composite=comp)
}

astro <- evaluate(get_candidates("Astro", c("Micro","Oligo")), "Astrocyte")
micro <- evaluate(get_candidates("Micro", c("Astro","Oligo")), "Microglia")
oligo <- evaluate(get_candidates("Oligo", c("Astro","Micro")), "Oligodendrocyte")

# Assemble covariate table
covtab <- data.frame(sample=sample_order, cohort=c(rep("2024",6),rep("2026",3)))
if (!is.null(astro)) covtab$astro <- round(astro$composite,4)
if (!is.null(micro)) covtab$micro <- round(micro$composite,4)
if (!is.null(oligo)) covtab$oligo <- round(oligo$composite,4)
cat("\n══ Covariate table ══\n"); print(covtab, row.names=FALSE)
write.table(covtab, "scMD_covariates.tsv", sep="\t", quote=FALSE, row.names=FALSE)
