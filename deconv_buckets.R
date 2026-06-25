suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
suppressPackageStartupMessages(library(SummarizedExperiment))

# Load HiBED libraries and extract layer probe IDs
data(HiBED_Libraries)
l1  <- rownames(HiBED_Libraries$Library_Layer1)
l2b <- rownames(HiBED_Libraries$Library_Layer2B)

# Annotation packages
anno_hg19 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_v2   <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2_base <- sub("_.*$", "", rownames(anno_v2))

chain <- import.chain("/projects2/DiffMeth/hg19ToHg38.over.chain")

# Build hybrid hg38 coordinate table for a set of probe IDs
# Uses EPICv2 native hg38 coords where available, liftover for remainder
get_hg38_coords <- function(probes, anno_hg19, anno_v2, anno_v2_base, chain) {
    # EPICv2 native coords for probes that match
    in_v2  <- probes[probes %in% anno_v2_base]
    v2_idx <- match(in_v2, anno_v2_base)
    coords_v2 <- data.frame(
        probe_id = in_v2,
        chr      = as.character(anno_v2$chr[v2_idx]),
        pos      = anno_v2$pos[v2_idx],
        stringsAsFactors = FALSE
    )

    # Liftover for probes not in EPICv2
    not_v2 <- probes[!probes %in% anno_v2_base]
    not_v2_in_hg19 <- not_v2[not_v2 %in% rownames(anno_hg19)]
    coords_lift <- NULL
    if (length(not_v2_in_hg19) > 0) {
        gr <- GRanges(
            seqnames  = anno_hg19[not_v2_in_hg19, "chr"],
            ranges    = IRanges(start = anno_hg19[not_v2_in_hg19, "pos"],
                                end   = anno_hg19[not_v2_in_hg19, "pos"]),
            probe_id  = not_v2_in_hg19
        )
        gr_hg38 <- unlist(liftOver(gr, chain))
        coords_lift <- data.frame(
            probe_id = gr_hg38$probe_id,
            chr      = as.character(seqnames(gr_hg38)),
            pos      = start(gr_hg38),
            stringsAsFactors = FALSE
        )
    }

    coords <- rbind(coords_v2, coords_lift)
    coords$coord <- paste(coords$chr, coords$pos, sep=":")
    coords
}

# Get hg38 coords for Layer 1 and Layer 2B
l1_coords  <- get_hg38_coords(l1,  anno_hg19, anno_v2, anno_v2_base, chain)
l2b_coords <- get_hg38_coords(l2b, anno_hg19, anno_v2, anno_v2_base, chain)

cat("Layer 1 probes with hg38 coords:", nrow(l1_coords), "\n")
cat("Layer 2B probes with hg38 coords:", nrow(l2b_coords), "\n")

# Sample file paths
bucket1_files <- c(
    CA192 = "CA192/CA192_GRCh38_CpG.bed",
    CA346 = "CA346/CA346_GRCh38_CpG.bed",
    CB239 = "CB239/CB239_GRCh38_CpG.bed",
    CC249 = "CC249/CC249_GRCh38_CpG.bed",
    CE167 = "CE167/CE167_GRCh38_CpG.bed",
    CE234 = "CE234/CE234_GRCh38_CpG.bed"
)

bucket2_files <- c(
    LG30 = "LG30/LG30_GRCh38_CpG.bed",
    LG31 = "LG31/LG31_GRCh38_CpG.bed",
    LG52 = "LG52/LG52_GRCh38_CpG.bed"
)

# Pool samples within a bucket: coverage-weighted mean beta per marker
pool_bucket <- function(bed_files, probe_coords) {
    all_betas <- lapply(names(bed_files), function(sname) {
        bed <- fread(bed_files[[sname]],
                     col.names=c("chr","start","end","beta","coverage"))
        bed$coord <- paste(bed$chr, bed$start, sep=":")
        merged <- merge(
            probe_coords[, .(coord, probe_id)],
            bed[, .(coord, beta, coverage)],
            by="coord", all.x=TRUE
        )
        merged$sample <- sname
        merged
    })
    combined <- rbindlist(all_betas)
    pooled <- combined[!is.na(beta), .(
        beta_pooled = weighted.mean(beta, coverage),
        n_samples   = .N
    ), by=.(probe_id)]
    pooled
}

b1_l1  <- pool_bucket(bucket1_files, as.data.table(l1_coords))
b2_l1  <- pool_bucket(bucket2_files, as.data.table(l1_coords))
b1_l2b <- pool_bucket(bucket1_files, as.data.table(l2b_coords))
b2_l2b <- pool_bucket(bucket2_files, as.data.table(l2b_coords))

cat("\nBucket 1 Layer 1 markers covered:", nrow(b1_l1), "\n")
cat("Bucket 2 Layer 1 markers covered:", nrow(b2_l1), "\n")
cat("Bucket 1 Layer 2B markers covered:", nrow(b1_l2b), "\n")
cat("Bucket 2 Layer 2B markers covered:", nrow(b2_l2b), "\n")

# Build beta matrix: rows=probes, cols=buckets
make_beta_matrix <- function(b1, b2, all_probes) {
    m <- matrix(NA, nrow=length(all_probes), ncol=2,
                dimnames=list(all_probes, c("Bucket1_2024","Bucket2_2026")))
    m[b1$probe_id, "Bucket1_2024"] <- b1$beta_pooled
    m[b2$probe_id, "Bucket2_2026"] <- b2$beta_pooled
    m
}

beta_l1  <- make_beta_matrix(b1_l1,  b2_l1,  l1_coords$probe_id)
beta_l2b <- make_beta_matrix(b1_l2b, b2_l2b, l2b_coords$probe_id)

cat("\nLayer 1 beta matrix dimensions:", dim(beta_l1), "\n")
cat("Non-NA values:", sum(!is.na(beta_l1)), "\n")
cat("\nLayer 2B beta matrix dimensions:", dim(beta_l2b), "\n")
cat("Non-NA values:", sum(!is.na(beta_l2b)), "\n")

# Run Layer 1 deconvolution
Library_Layer1 <- as.data.frame(assay(HiBED_Libraries$Library_Layer1, "counts"))
common_l1 <- intersect(rownames(beta_l1)[rowSums(!is.na(beta_l1)) > 0],
                       rownames(Library_Layer1))
cat("\nLayer 1 common probes for deconvolution:", length(common_l1), "\n")

result_l1 <- projectCellType_CP(
    beta_l1[common_l1, , drop=FALSE],
    as.matrix(Library_Layer1[common_l1, ]),
    lessThanOne=TRUE
)
cat("\nLayer 1 deconvolution results:\n")
print(result_l1)

# Run Layer 2B deconvolution
Library_Layer2B <- as.data.frame(assay(HiBED_Libraries$Library_Layer2B, "counts"))
common_l2b <- intersect(rownames(beta_l2b)[rowSums(!is.na(beta_l2b)) > 0],
                        rownames(Library_Layer2B))
cat("\nLayer 2B common probes for deconvolution:", length(common_l2b), "\n")

result_l2b <- projectCellType_CP(
    beta_l2b[common_l2b, , drop=FALSE],
    as.matrix(Library_Layer2B[common_l2b, ]),
    lessThanOne=TRUE
)
cat("\nLayer 2B deconvolution results:\n")
print(result_l2b)
