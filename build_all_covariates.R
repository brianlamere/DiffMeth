suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SummarizedExperiment))

setwd("/projects2/DiffMeth")
data(HiBED_Libraries)

anno_hg19    <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_v2      <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2_base <- sub("_.*$", "", rownames(anno_v2))
chain        <- import.chain("/projects2/DiffMeth/hg19ToHg38.over.chain")

get_hg38_coords <- function(probes, anno_hg19, anno_v2, anno_v2_base, chain) {
    in_v2  <- probes[probes %in% anno_v2_base]
    v2_idx <- match(in_v2, anno_v2_base)
    coords_v2 <- data.frame(probe_id=in_v2,
        chr=as.character(anno_v2$chr[v2_idx]), pos=anno_v2$pos[v2_idx],
        stringsAsFactors=FALSE)
    not_v2 <- probes[!probes %in% anno_v2_base]
    not_v2_in_hg19 <- not_v2[not_v2 %in% rownames(anno_hg19)]
    coords_lift <- NULL
    if (length(not_v2_in_hg19) > 0) {
        gr <- GRanges(seqnames=anno_hg19[not_v2_in_hg19,"chr"],
            ranges=IRanges(start=anno_hg19[not_v2_in_hg19,"pos"],
                           end=anno_hg19[not_v2_in_hg19,"pos"]),
            probe_id=not_v2_in_hg19)
        gr_hg38 <- unlist(liftOver(gr, chain))
        coords_lift <- data.frame(probe_id=gr_hg38$probe_id,
            chr=as.character(seqnames(gr_hg38)), pos=start(gr_hg38),
            stringsAsFactors=FALSE)
    }
    coords <- rbind(coords_v2, coords_lift)
    coords$coord <- paste(coords$chr, coords$pos, sep=":")
    coords
}

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
sample_beds <- setNames(
    file.path(sample_order, paste0(sample_order, "_GRCh38_CpG.bed")), sample_order)
bed_list <- lapply(sample_beds, function(f) {
    b <- fread(f, col.names=c("chr","start","end","beta","coverage"))
    b$coord <- paste(b$chr, b$start, sep=":"); setkey(b, coord); b
})

ref2b <- as.data.frame(assay(HiBED_Libraries$Library_Layer2B, "counts"))
l2b   <- rownames(ref2b)

# For one cell type: select usable markers, return coverage-weighted composite per sample
build_covariate <- function(target_col, label) {
    other <- setdiff(colnames(ref2b), target_col)
    delta <- ref2b[l2b, target_col] - rowMeans(ref2b[l2b, other, drop=FALSE])
    specific <- l2b[abs(delta) > 0.3]
    coords <- get_hg38_coords(specific, anno_hg19, anno_v2, anno_v2_base, chain)

    cov  <- sapply(sample_order, function(s) bed_list[[s]][.(coords$coord), coverage])
    beta <- sapply(sample_order, function(s) bed_list[[s]][.(coords$coord), beta])
    rownames(cov) <- coords$probe_id; rownames(beta) <- coords$probe_id

    # usable = covered in all 9, sd>0.1, all depths >=4
    usable <- character(0)
    for (i in seq_len(nrow(cov))) {
        if (all(!is.na(cov[i,])) && sd(beta[i,]) > 0.1 && all(cov[i,] >= 4))
            usable <- c(usable, rownames(cov)[i])
    }
    cat(sprintf("\n%s: %d usable markers: %s\n", label, length(usable),
                paste(usable, collapse=", ")))
    if (length(usable) == 0) return(rep(NA, length(sample_order)))

    # coverage-weighted composite across usable markers, per sample
    comp <- sapply(sample_order, function(s) {
        b <- beta[usable, s]; w <- cov[usable, s]
        sum(b*w) / sum(w)
    })
    comp
}

cov_oligo  <- build_covariate("Oligodendrocyte", "Oligodendrocyte")
cov_astro  <- build_covariate("Astrocyte",       "Astrocyte")
cov_micro  <- build_covariate("Microglial",      "Microglial")
cov_neuron <- build_covariate("Neuronal",        "Neuronal")

covariates <- data.frame(
    sample = sample_order,
    cohort = c(rep("2024",6), rep("2026",3)),
    oligo  = round(cov_oligo, 4),
    astro  = round(cov_astro, 4),
    micro  = round(cov_micro, 4),
    neuron = round(cov_neuron, 4)
)

cat("\n══ Per-sample composite covariates ══\n")
print(covariates, row.names=FALSE)

cat("\n══ Per-cohort means (which cell types differ most between buckets) ══\n")
for (ct in c("oligo","astro","micro","neuron")) {
    m24 <- mean(covariates[covariates$cohort=="2024", ct], na.rm=TRUE)
    m26 <- mean(covariates[covariates$cohort=="2026", ct], na.rm=TRUE)
    cat(sprintf("  %-8s 2024=%.3f  2026=%.3f  gap=%.3f\n", ct, m24, m26, abs(m24-m26)))
}

write.table(covariates, "composition_covariates.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWritten: composition_covariates.tsv\n")
