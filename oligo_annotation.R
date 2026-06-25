suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(SummarizedExperiment))

setwd("/projects2/DiffMeth")

# ══ 1. Load the three result sets ════════════════════════════════════════════
sig_unadj <- fread("sig_unadj.tsv")
sig_covA  <- fread("sig_adj_cg16515546.tsv")
sig_covB  <- fread("sig_adj_composite.tsv")

key <- function(dt) paste(dt$chr, dt$start, dt$end)
ku <- key(sig_unadj); kA <- key(sig_covA); kB <- key(sig_covB)

# Classify unadjusted tiles: lost under BOTH adjustments = most clearly composition-linked
sig_unadj[, survived_A := ku %in% kA]
sig_unadj[, survived_B := ku %in% kB]
sig_unadj[, status := fifelse(survived_A & survived_B, "survived_both",
                       fifelse(!survived_A & !survived_B, "lost_both", "mixed"))]

cat("══ Unadjusted tile fate after adjustment ══\n")
print(table(sig_unadj$status))

# ══ 2. Build oligodendrocyte-specific position set from HiBED Layer 2B ════════
data(HiBED_Libraries)
ref2b <- as.data.frame(assay(HiBED_Libraries$Library_Layer2B, "counts"))
l2b   <- rownames(ref2b)

oligo_col   <- "Oligodendrocyte"
other_types <- setdiff(colnames(ref2b), oligo_col)
oligo_delta <- ref2b[[oligo_col]] - rowMeans(ref2b[, other_types])
# Oligodendrocyte-specific markers: strongly distinct in either direction
oligo_probes <- l2b[abs(oligo_delta) > 0.3]
cat("\nOligodendrocyte-specific Layer 2B markers:", length(oligo_probes), "\n")

# ══ 3. Coordinate conversion ═════════════════════════════════════════════════
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
    rbind(coords_v2, coords_lift)
}

oligo_coords <- get_hg38_coords(oligo_probes, anno_hg19, anno_v2, anno_v2_base, chain)
cat("Oligo markers with hg38 coords:", nrow(oligo_coords), "\n")

# ══ 4. Overlap each tile set with oligodendrocyte positions ══════════════════
oligo_gr <- GRanges(seqnames=oligo_coords$chr,
                    ranges=IRanges(start=oligo_coords$pos, end=oligo_coords$pos))

tile_overlap <- function(dt) {
    gr <- GRanges(seqnames=dt$chr, ranges=IRanges(start=dt$start, end=dt$end))
    countOverlaps(gr, oligo_gr) > 0
}

sig_unadj[, oligo_overlap := tile_overlap(sig_unadj)]

# ══ 5. The key test: are LOST tiles preferentially oligodendrocyte-associated? ═
cat("\n══ Oligodendrocyte-marker overlap by tile fate ══\n")
tab <- table(status = sig_unadj$status, oligo = sig_unadj$oligo_overlap)
print(tab)

cat("\nProportion overlapping oligo markers, by fate:\n")
for (st in unique(sig_unadj$status)) {
    sub <- sig_unadj[status == st]
    cat(sprintf("  %-14s %d/%d (%.1f%%)\n", st,
        sum(sub$oligo_overlap), nrow(sub),
        100*mean(sub$oligo_overlap)))
}

# Fisher test: lost_both vs survived_both, oligo overlap enrichment
lost  <- sig_unadj[status=="lost_both"]
surv  <- sig_unadj[status=="survived_both"]
if (nrow(lost) > 0 && nrow(surv) > 0) {
    ft <- fisher.test(matrix(c(
        sum(lost$oligo_overlap),  sum(!lost$oligo_overlap),
        sum(surv$oligo_overlap),  sum(!surv$oligo_overlap)), nrow=2))
    cat(sprintf("\nFisher test (lost_both vs survived_both oligo enrichment):\n"))
    cat(sprintf("  odds ratio = %.2f, p = %.4f\n", ft$estimate, ft$p.value))
}

# ══ 6. Save annotated table ══════════════════════════════════════════════════
write.table(sig_unadj, "sig_unadj_oligo_annotated.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWritten: sig_unadj_oligo_annotated.tsv\n")
