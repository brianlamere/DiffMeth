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

# Pre-load all sample BEDs once (keyed by coord) for speed
bed_list <- lapply(sample_beds, function(f) {
    b <- fread(f, col.names=c("chr","start","end","beta","coverage"))
    b$coord <- paste(b$chr, b$start, sep=":")
    setkey(b, coord)
    b
})

# For a given marker set + a "target" cell type column, find usable markers
scan_celltype <- function(probes, ref, target_col, label) {
    other <- setdiff(colnames(ref), target_col)
    delta <- ref[probes, target_col] - rowMeans(ref[probes, other, drop=FALSE])
    specific <- probes[abs(delta) > 0.3]
    if (length(specific) == 0) {
        cat(sprintf("%-14s no markers with |delta|>0.3\n", label)); return(invisible())
    }
    coords <- get_hg38_coords(specific, anno_hg19, anno_v2, anno_v2_base, chain)
    # coverage + beta across samples
    cov <- sapply(sample_order, function(s) {
        bed_list[[s]][.(coords$coord), coverage]
    })
    beta <- sapply(sample_order, function(s) {
        bed_list[[s]][.(coords$coord), beta]
    })
    rownames(cov) <- coords$probe_id; rownames(beta) <- coords$probe_id
    universal <- which(rowSums(!is.na(cov)) == length(sample_order))
    n_univ <- length(universal)
    # of universal, how many have dynamic range (sd > 0.1) and no single-read extremes
    usable <- 0
    if (n_univ > 0) {
        for (i in universal) {
            sd_ok    <- sd(beta[i, ], na.rm=TRUE) > 0.1
            depth_ok <- all(cov[i, ] >= 4)   # no value resting on <4 reads
            if (sd_ok && depth_ok) usable <- usable + 1
        }
    }
    cat(sprintf("%-14s specific:%3d  universal:%2d  usable(range+depth):%2d\n",
                label, length(specific), n_univ, usable))
}

# Layer 1: Neuronal / Glial / Endothelial-Stromal
ref1 <- as.data.frame(assay(HiBED_Libraries$Library_Layer1, "counts"))
l1 <- rownames(ref1)
cat("══ Layer 1 ══\n")
for (ct in colnames(ref1)) scan_celltype(l1, ref1, ct, ct)

# Layer 2B: Astrocyte / Microglial / Oligodendrocyte / Neuronal
ref2b <- as.data.frame(assay(HiBED_Libraries$Library_Layer2B, "counts"))
l2b <- rownames(ref2b)
cat("\n══ Layer 2B ══\n")
for (ct in colnames(ref2b)) scan_celltype(l2b, ref2b, ct, ct)

# Layer 2C: GABA / GLU neuronal subtypes
ref2c <- as.data.frame(assay(HiBED_Libraries$Library_Layer2C, "counts"))
l2c <- rownames(ref2c)
cat("\n══ Layer 2C ══\n")
for (ct in colnames(ref2c)) scan_celltype(l2c, ref2c, ct, ct)
