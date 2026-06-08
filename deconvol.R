suppressPackageStartupMessages(library(CETYGO))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2manifest))
suppressPackageStartupMessages(library(data.table))

data(HiBED_Libraries)

# Get all unique probe IDs across all IDOL sub-models
all_probes <- unique(unlist(lapply(modelBrainCoef$IDOL, rownames)))

# Get hg19 coordinates for these probes
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
marker_anno <- anno[all_probes[all_probes %in% rownames(anno)], 
                    c("chr","pos")]
marker_anno <- marker_anno[!is.na(marker_anno$chr), ]

# Now liftover hg19 -> hg38
chain <- import.chain("/projects2/DiffMeth/hg19ToHg38.over.chain")

marker_gr <- GRanges(
    seqnames = marker_anno$chr,
    ranges = IRanges(start=marker_anno$pos, end=marker_anno$pos),
    probe_id = rownames(marker_anno)
)

marker_hg38 <- unlist(liftOver(marker_gr, chain))

# Now check against your RRBS coverage - use CA192 as representative
ca192 <- read.table("CA192/CA192_GRCh38_CpG.bed",
                    col.names=c("chr","start","end","beta","coverage"))
ca192$coord <- paste(ca192$chr, ca192$start, sep=":")

marker_coords_hg38 <- paste(as.character(seqnames(marker_hg38)),
                             start(marker_hg38), sep=":")

overlap <- sum(marker_coords_hg38 %in% ca192$coord)
cat("Total unique marker CpGs:", length(all_probes), "\n")
cat("With hg19 coordinates:", nrow(marker_anno), "\n") 
cat("Survived liftover to hg38:", length(marker_hg38), "\n")
cat("Covered in your RRBS:", overlap, "\n")
cat("Overlap %:", round(overlap/length(marker_hg38)*100, 1), "%\n")

l1 <- rownames(HiBED_Libraries$Library_Layer1)
l2a <- rownames(HiBED_Libraries$Library_Layer2A)
l2b <- rownames(HiBED_Libraries$Library_Layer2B)
l2c <- rownames(HiBED_Libraries$Library_Layer2C)

cat("Layer 1 markers:", length(l1), "\n")
cat("Layer 2A markers:", length(l2a), "\n")
cat("Layer 2B markers:", length(l2b), "\n")
cat("Layer 2C markers:", length(l2c), "\n")

# Look at which specific Layer 1 markers you have
l1_anno <- anno[l1[l1 %in% rownames(anno)], c("chr","pos")]
l1_gr <- GRanges(seqnames=l1_anno$chr,
                 ranges=IRanges(start=l1_anno$pos, end=l1_anno$pos),
                 probe_id=rownames(l1_anno))
l1_hg38 <- unlist(liftOver(l1_gr, chain))
l1_coords <- paste(as.character(seqnames(l1_hg38)),
                   start(l1_hg38), sep=":")

# Extract probe IDs from metadata column
l1_coords_hg38 <- paste(as.character(seqnames(l1_hg38)),
                         start(l1_hg38), sep=":")

# Look at which specific Layer 1 markers you have
l2b_anno <- anno[l2b[l2b %in% rownames(anno)], c("chr","pos")]
l2b_gr <- GRanges(seqnames=l2b_anno$chr,
                 ranges=IRanges(start=l2b_anno$pos, end=l2b_anno$pos),
                 probe_id=rownames(l2b_anno))
l2b_hg38 <- unlist(liftOver(l2b_gr, chain))
l2b_coords <- paste(as.character(seqnames(l2b_hg38)),
                   start(l2b_hg38), sep=":")

# Extract probe IDs from metadata column
l2b_coords_hg38 <- paste(as.character(seqnames(l2b_hg38)),
                         start(l2b_hg38), sep=":")

bed_file <- list("CA192/CA192_GRCh38_CpG.bed")
ca192_rrbs <- rbindlist(lapply(bed_file, function(f) {
    fread(f, col.names=c("chr","start","end","beta","coverage"))
}))
rrbs_coords <- unique(paste(ca192_rrbs$chr, ca192_rrbs$start, sep=":"))
cat("Total unique CpGs across sample:", length(rrbs_coords), "\n")

# Match against your RRBS coords
covered_idx <- l1_coords_hg38 %in% rrbs_coords
recovered_probes <- l1_hg38$probe_id[covered_idx]
cat("Recovered Layer 1 probes:\n")
print(recovered_probes)
cat("Count:", length(recovered_probes), "\n")

assay(HiBED_Libraries$Library_Layer1)[recovered_probes, ]

# Match against your RRBS coords
covered_idx_2b <- l2b_coords_hg38 %in% rrbs_coords
recovered_probes_2b <- l2b_hg38$probe_id[covered_idx]
cat("Recovered Layer 2b probes:\n")
print(recovered_probes_2b)
cat("Count:", length(recovered_probes_2b), "\n")

assay(HiBED_Libraries$Library_Layer2B)[recovered_probes_2b, ]

---------------

# Which ones did you recover?
recovered_probes <- names(l1_hg38)[l1_coords %in% rrbs_coords]
cat("Recovered Layer 1 probes:\n")
print(recovered_probes)

# Look at their methylation values in the reference
assay(HiBED_Libraries$Library_Layer1)[recovered_probes, ]


# Read each sample's CpG bed and extract beta for those coordinates
sample_beds <- list(
    CA192 = "CA192/CA192_GRCh38_CpG.bed",
    CA346 = "CA346/CA346_GRCh38_CpG.bed",
    CB239 = "CB239/CB239_GRCh38_CpG.bed",
    CC249 = "CC249/CC249_GRCh38_CpG.bed",
    CE167 = "CE167/CE167_GRCh38_CpG.bed",
    CE234 = "CE234/CE234_GRCh38_CpG.bed",
    LG30  = "LG30/LG30_GRCh38_CpG.bed",
    LG31  = "LG31/LG31_GRCh38_CpG.bed",
    LG52  = "LG52/LG52_GRCh38_CpG.bed"
)


# What coordinate did cg15322932 lift over to?
idx <- which(l1_hg38$probe_id == "cg15322932")
cat("hg38 coordinate:", as.character(seqnames(l1_hg38[idx])),
    start(l1_hg38[idx]), "\n")

# What does CA192 have near that position?
ca192 <- fread("CA192/CA192_GRCh38_CpG.bed",
               col.names=c("chr","start","end","beta","coverage"))
ca192$coord <- paste(ca192$chr, ca192$start, sep=":")

# Check the exact coordinate
target_coord <- paste(as.character(seqnames(l1_hg38[idx])),
                      start(l1_hg38[idx]), sep=":")
cat("Looking for:", target_coord, "\n")
cat("Found in CA192:", target_coord %in% ca192$coord, "\n")

# Look at nearby positions in case of off-by-one
target_pos <- start(l1_hg38[idx])
ca192[abs(ca192$start - target_pos) < 5 &
      ca192$chr == as.character(seqnames(l1_hg38[idx])), ]
