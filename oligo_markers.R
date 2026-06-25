suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SummarizedExperiment))

setwd("/projects2/DiffMeth")

data(HiBED_Libraries)
l2b <- rownames(HiBED_Libraries$Library_Layer2B)

# Reference methylation values for Layer 2B (cell types in columns)
ref2b <- as.data.frame(assay(HiBED_Libraries$Library_Layer2B, "counts"))
cat("Layer 2B cell types:", paste(colnames(ref2b), collapse=", "), "\n\n")

# ── Coordinate conversion (same hybrid approach) ──────────────────────────────
anno_hg19    <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_v2      <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2_base <- sub("_.*$", "", rownames(anno_v2))
chain        <- import.chain("/projects2/DiffMeth/hg19ToHg38.over.chain")

get_hg38_coords <- function(probes, anno_hg19, anno_v2, anno_v2_base, chain) {
    in_v2    <- probes[probes %in% anno_v2_base]
    v2_idx   <- match(in_v2, anno_v2_base)
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

l2b_coords <- get_hg38_coords(l2b, anno_hg19, anno_v2, anno_v2_base, chain)

# ── Sample BEDs ───────────────────────────────────────────────────────────────
sample_beds <- c(
    CA192="CA192/CA192_GRCh38_CpG.bed", CA346="CA346/CA346_GRCh38_CpG.bed",
    CB239="CB239/CB239_GRCh38_CpG.bed", CC249="CC249/CC249_GRCh38_CpG.bed",
    CE167="CE167/CE167_GRCh38_CpG.bed", CE234="CE234/CE234_GRCh38_CpG.bed",
    LG30="LG30/LG30_GRCh38_CpG.bed", LG31="LG31/LG31_GRCh38_CpG.bed",
    LG52="LG52/LG52_GRCh38_CpG.bed")

# Build coverage + beta matrices
coord_set <- as.data.table(l2b_coords)
cov_mat  <- matrix(NA, nrow=nrow(coord_set), ncol=length(sample_beds),
                   dimnames=list(coord_set$probe_id, names(sample_beds)))
beta_mat <- cov_mat
for (s in names(sample_beds)) {
    bed <- fread(sample_beds[[s]], col.names=c("chr","start","end","beta","coverage"))
    bed$coord <- paste(bed$chr, bed$start, sep=":")
    m <- match(coord_set$coord, bed$coord)
    cov_mat[, s]  <- bed$coverage[m]
    beta_mat[, s] <- bed$beta[m]
}

# Universally covered markers
universal <- rownames(cov_mat)[rowSums(!is.na(cov_mat)) == ncol(cov_mat)]
cat("Universally-covered Layer 2B markers:", length(universal), "\n\n")

# ── Oligodendrocyte specificity in the reference ─────────────────────────────
# How distinct is Oligodendrocyte from the other 3 cell types at each marker?
oligo_col   <- "Oligodendrocyte"
other_types <- setdiff(colnames(ref2b), oligo_col)

oligo_spec <- data.frame(
    probe_id    = universal,
    oligo_ref   = ref2b[universal, oligo_col],
    others_mean = rowMeans(ref2b[universal, other_types]),
    stringsAsFactors = FALSE
)
oligo_spec$oligo_delta <- oligo_spec$oligo_ref - oligo_spec$others_mean

# Per-sample observed beta spread (dynamic range) in YOUR data
oligo_spec$obs_min  <- apply(beta_mat[universal, ], 1, min)
oligo_spec$obs_max  <- apply(beta_mat[universal, ], 1, max)
oligo_spec$obs_range <- oligo_spec$obs_max - oligo_spec$obs_min
oligo_spec$obs_sd    <- apply(beta_mat[universal, ], 1, sd)

# Within-cohort vs between-cohort: do the betas track the 2024/2026 split?
b2024 <- c("CA192","CA346","CB239","CC249","CE167","CE234")
b2026 <- c("LG30","LG31","LG52")
oligo_spec$mean_2024 <- rowMeans(beta_mat[universal, b2024])
oligo_spec$mean_2026 <- rowMeans(beta_mat[universal, b2026])
oligo_spec$cohort_gap <- abs(oligo_spec$mean_2024 - oligo_spec$mean_2026)

# Sort by oligodendrocyte specificity (strongest first)
oligo_spec <- oligo_spec[order(-abs(oligo_spec$oligo_delta)), ]

cat("Layer 2B universally-covered markers, ranked by oligodendrocyte specificity:\n")
cat("(oligo_delta = how distinct Oligo is from other cell types in reference)\n")
cat("(obs_range/obs_sd = dynamic range in YOUR samples)\n")
cat("(cohort_gap = how much betas separate by cohort - LOW is better, avoids collinearity)\n\n")
print(round(oligo_spec[, c("oligo_ref","others_mean","oligo_delta",
                            "obs_range","obs_sd","mean_2024","mean_2026",
                            "cohort_gap")], 3))

# Show per-sample betas for the top oligodendrocyte-specific markers
top_oligo <- head(oligo_spec$probe_id[abs(oligo_spec$oligo_delta) > 0.3], 6)
cat("\nPer-sample betas for strongest oligodendrocyte-specific markers:\n")
print(round(beta_mat[top_oligo, , drop=FALSE], 3))
