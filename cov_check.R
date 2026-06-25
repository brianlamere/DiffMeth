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

markers <- c("cg16515546", "cg12967384")

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

mc <- get_hg38_coords(markers, anno_hg19, anno_v2, anno_v2_base, chain)
print(mc)

sample_beds <- c(
    CA192="CA192/CA192_GRCh38_CpG.bed", CA346="CA346/CA346_GRCh38_CpG.bed",
    CB239="CB239/CB239_GRCh38_CpG.bed", CC249="CC249/CC249_GRCh38_CpG.bed",
    CE167="CE167/CE167_GRCh38_CpG.bed", CE234="CE234/CE234_GRCh38_CpG.bed",
    LG30="LG30/LG30_GRCh38_CpG.bed", LG31="LG31/LG31_GRCh38_CpG.bed",
    LG52="LG52/LG52_GRCh38_CpG.bed")

cat("\nCoverage depth (beta) per sample:\n\n")
for (i in seq_len(nrow(mc))) {
    cat("══", mc$probe_id[i], "(", mc$coord[i], ") ══\n")
    for (s in names(sample_beds)) {
        bed <- fread(sample_beds[[s]], col.names=c("chr","start","end","beta","coverage"))
        bed$coord <- paste(bed$chr, bed$start, sep=":")
        row <- bed[coord == mc$coord[i]]
        if (nrow(row) > 0) {
            cat(sprintf("  %-6s depth=%-4d beta=%.3f\n", s, row$coverage[1], row$beta[1]))
        } else {
            cat(sprintf("  %-6s NOT COVERED\n", s))
        }
    }
    cat("\n")
}

