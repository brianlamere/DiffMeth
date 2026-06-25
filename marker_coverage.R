suppressPackageStartupMessages(library(HiBED))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SummarizedExperiment))

setwd("/projects2/DiffMeth")

# ── HiBED marker probe IDs ────────────────────────────────────────────────────
data(HiBED_Libraries)
l1  <- rownames(HiBED_Libraries$Library_Layer1)
l2b <- rownames(HiBED_Libraries$Library_Layer2B)

# ── Hybrid hg38 coordinate conversion ────────────────────────────────────────
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
    not_v2         <- probes[!probes %in% anno_v2_base]
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

# ── Sample BED files ──────────────────────────────────────────────────────────
sample_beds <- c(
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

# ── Build coverage matrix: rows=markers, cols=samples ────────────────────────
# Value = coverage depth at that marker in that sample (NA if not covered)
build_coverage_matrix <- function(probe_coords, sample_beds) {
    coord_set <- as.data.table(probe_coords)
    mat <- matrix(NA, nrow=nrow(coord_set), ncol=length(sample_beds),
                  dimnames=list(coord_set$probe_id, names(sample_beds)))
    beta_mat <- mat
    for (sname in names(sample_beds)) {
        bed <- fread(sample_beds[[sname]],
                     col.names=c("chr","start","end","beta","coverage"))
        bed$coord <- paste(bed$chr, bed$start, sep=":")
        m <- match(coord_set$coord, bed$coord)
        mat[, sname]      <- bed$coverage[m]
        beta_mat[, sname] <- bed$beta[m]
    }
    list(coverage=mat, beta=beta_mat)
}

l1_cov  <- build_coverage_matrix(l1_coords,  sample_beds)
l2b_cov <- build_coverage_matrix(l2b_coords, sample_beds)

# ── Report: markers covered in all 9 samples ──────────────────────────────────
report_layer <- function(cov, beta, label) {
    n_covered  <- rowSums(!is.na(cov))
    all_nine   <- which(n_covered == ncol(cov))

    cat("\n══ ", label, " ══\n", sep="")
    cat("Total markers:", nrow(cov), "\n")
    cat("Distribution of how many samples cover each marker:\n")
    print(table(n_covered))
    cat("\nMarkers covered in ALL 9 samples:", length(all_nine), "\n")

    if (length(all_nine) > 0) {
        cat("\nUniversally-covered markers (coverage depth per sample):\n")
        print(cov[all_nine, , drop=FALSE])
        cat("\nUniversally-covered markers (methylation beta per sample):\n")
        print(round(beta[all_nine, , drop=FALSE], 3))
    }
    invisible(all_nine)
}

report_layer(l1_cov$coverage,  l1_cov$beta,  "HiBED Layer 1 (Neuronal / Glial / Endothelial-Stromal)")
report_layer(l2b_cov$coverage, l2b_cov$beta, "HiBED Layer 2B (Astrocyte / Microglial / Oligodendrocyte / Neuronal)")
