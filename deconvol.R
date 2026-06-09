suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
suppressPackageStartupMessages(library(SummarizedExperiment))

setwd("/projects2/DiffMeth")

# ── 1. HiBED layer probe IDs ──────────────────────────────────────────────────
data(HiBED_Libraries)
l1  <- rownames(HiBED_Libraries$Library_Layer1)
l2b <- rownames(HiBED_Libraries$Library_Layer2B)

# ── 2. Annotation: hybrid hg38 coordinates ───────────────────────────────────
# EPICv2 has native hg38; strip suffix to match probe IDs
# Fall back to hg19 liftover for probes absent from EPICv2
anno_hg19    <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_v2      <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2_base <- sub("_.*$", "", rownames(anno_v2))
chain        <- import.chain("/projects2/DiffMeth/hg19ToHg38.over.chain")

get_hg38_coords <- function(probes, anno_hg19, anno_v2, anno_v2_base, chain) {
    in_v2    <- probes[probes %in% anno_v2_base]
    v2_idx   <- match(in_v2, anno_v2_base)
    coords_v2 <- data.frame(
        probe_id = in_v2,
        chr      = as.character(anno_v2$chr[v2_idx]),
        pos      = anno_v2$pos[v2_idx],
        stringsAsFactors = FALSE
    )
    not_v2        <- probes[!probes %in% anno_v2_base]
    not_v2_in_hg19 <- not_v2[not_v2 %in% rownames(anno_hg19)]
    coords_lift <- NULL
    if (length(not_v2_in_hg19) > 0) {
        gr <- GRanges(
            seqnames = anno_hg19[not_v2_in_hg19, "chr"],
            ranges   = IRanges(start = anno_hg19[not_v2_in_hg19, "pos"],
                               end   = anno_hg19[not_v2_in_hg19, "pos"]),
            probe_id = not_v2_in_hg19
        )
        gr_hg38 <- unlist(liftOver(gr, chain))
        coords_lift <- data.frame(
            probe_id = gr_hg38$probe_id,
            chr      = as.character(seqnames(gr_hg38)),
            pos      = start(gr_hg38),
            stringsAsFactors = FALSE
        )
    }
    coords       <- rbind(coords_v2, coords_lift)
    coords$coord <- paste(coords$chr, coords$pos, sep=":")
    coords
}

l1_coords  <- get_hg38_coords(l1,  anno_hg19, anno_v2, anno_v2_base, chain)
l2b_coords <- get_hg38_coords(l2b, anno_hg19, anno_v2, anno_v2_base, chain)

cat("Layer 1 probes with hg38 coords:", nrow(l1_coords), "\n")
cat("Layer 2B probes with hg38 coords:", nrow(l2b_coords), "\n")

# ── 3. Sample file paths ──────────────────────────────────────────────────────
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

# ── 4. Pool betas within bucket (coverage-weighted mean) ─────────────────────
pool_bucket <- function(bed_files, probe_coords) {
    all_betas <- lapply(names(bed_files), function(sname) {
        bed <- fread(bed_files[[sname]],
                     col.names=c("chr","start","end","beta","coverage"))
        bed$coord <- paste(bed$chr, bed$start, sep=":")
        merged <- merge(
            as.data.table(probe_coords)[, .(coord, probe_id)],
            bed[, .(coord, beta, coverage)],
            by="coord", all.x=TRUE
        )
        merged$sample <- sname
        merged
    })
    combined <- rbindlist(all_betas)
    combined[!is.na(beta), .(
        beta_pooled = weighted.mean(beta, coverage),
        n_samples   = .N
    ), by=.(probe_id)]
}

b1_l1  <- pool_bucket(bucket1_files, l1_coords)
b2_l1  <- pool_bucket(bucket2_files, l1_coords)
b1_l2b <- pool_bucket(bucket1_files, l2b_coords)
b2_l2b <- pool_bucket(bucket2_files, l2b_coords)

# ── 5. Build beta matrices ────────────────────────────────────────────────────
make_beta_matrix <- function(b1, b2, all_probes) {
    m <- matrix(NA, nrow=length(all_probes), ncol=2,
                dimnames=list(all_probes, c("Bucket1_2024","Bucket2_2026")))
    m[b1$probe_id, "Bucket1_2024"] <- b1$beta_pooled
    m[b2$probe_id, "Bucket2_2026"] <- b2$beta_pooled
    m
}

beta_l1  <- make_beta_matrix(b1_l1,  b2_l1,  l1_coords$probe_id)
beta_l2b <- make_beta_matrix(b1_l2b, b2_l2b, l2b_coords$probe_id)

# ── 6. Identify shared markers (covered in both buckets) ─────────────────────
shared_l1  <- intersect(b1_l1$probe_id,  b2_l1$probe_id)
shared_l2b <- intersect(b1_l2b$probe_id, b2_l2b$probe_id)
cat("Layer 1 shared markers:", length(shared_l1), "\n")
cat("Layer 2B shared markers:", length(shared_l2b), "\n")

# ── 7. Deconvolve using shared markers only ───────────────────────────────────
Library_Layer1  <- as.data.frame(assay(HiBED_Libraries$Library_Layer1,  "counts"))
Library_Layer2B <- as.data.frame(assay(HiBED_Libraries$Library_Layer2B, "counts"))

result_l1 <- projectCellType_CP(
    beta_l1[shared_l1, , drop=FALSE],
    as.matrix(Library_Layer1[shared_l1, ]),
    lessThanOne = TRUE
)
cat("\nLayer 1 deconvolution (shared markers only):\n")
print(result_l1)

result_l2b <- projectCellType_CP(
    beta_l2b[shared_l2b, , drop=FALSE],
    as.matrix(Library_Layer2B[shared_l2b, ]),
    lessThanOne = TRUE
)
cat("\nLayer 2B deconvolution (shared markers only):\n")
print(result_l2b)

# ── 8. Summary for use as methylKit covariates ────────────────────────────────
cat("\n── Glial fractions for use as covariates ──\n")
cat("Bucket1_2024 Glial:", round(result_l1["Bucket1_2024","Glial"], 4), "\n")
cat("Bucket2_2026 Glial:", round(result_l1["Bucket2_2026","Glial"], 4), "\n")
cat("Bucket1_2024 Neuronal:", round(result_l1["Bucket1_2024","Neuronal"], 4), "\n")
cat("Bucket2_2026 Neuronal:", round(result_l1["Bucket2_2026","Neuronal"], 4), "\n")
