suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

setwd("/projects2/DiffMeth")

# ── 1. Load significant tiles ─────────────────────────────────────────────────
#sig_tiles <- read.table("sig_tiles_unadjusted.tsv", header=TRUE, sep="\t")
sig_tiles <- read.table("sig_tiles_unadjusted_q010.tsv", header=TRUE, sep="\t")
cat("Significant tiles loaded:", nrow(sig_tiles), "\n")

# ── 2. Load Loyfer 2023 methylation atlas ─────────────────────────────────────
header_line <- readLines(
    gzcon(file("/projects1/references/meth_atlas/full_atlas.csv.gz", "rb")), n=1)
col_names    <- strsplit(header_line, ",")[[1]]
col_names[1] <- "probe_id"

atlas <- fread("/projects1/references/meth_atlas/full_atlas.csv.gz",
               skip=1, header=FALSE)
setnames(atlas, col_names)
cat("Atlas probes loaded:", nrow(atlas), "\n")

# ── 3. Identify neuron-specific probes ────────────────────────────────────────
# Probes where Cortical_neurons methylation differs from mean of all other
# cell types by more than 0.3 in either direction
other_cols <- setdiff(names(atlas), c("probe_id", "Cortical_neurons"))
atlas[, mean_other  := rowMeans(.SD, na.rm=TRUE), .SDcols=other_cols]
atlas[, neuron_delta := Cortical_neurons - mean_other]

neuron_specific <- atlas[abs(neuron_delta) > 0.3, probe_id]
cat("Neuron-specific probes (|delta| > 0.3):", length(neuron_specific), "\n")

# ── 4. Hybrid hg38 coordinate conversion ─────────────────────────────────────
# EPICv2 for native hg38 coords; liftover for remainder
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

neuron_coords <- get_hg38_coords(neuron_specific,
                                  anno_hg19, anno_v2, anno_v2_base, chain)
cat("Neuron-specific probes with hg38 coords:", nrow(neuron_coords), "\n")

# ── 5. Intersect with significant tiles ───────────────────────────────────────
sig_gr <- GRanges(
    seqnames = sig_tiles$chr,
    ranges   = IRanges(start=sig_tiles$start, end=sig_tiles$end)
)

neuron_gr <- GRanges(
    seqnames = neuron_coords$chr,
    ranges   = IRanges(start=neuron_coords$pos, end=neuron_coords$pos)
)

hits <- findOverlaps(sig_gr, neuron_gr)
tiles_with_neuron_signal <- unique(queryHits(hits))

# ── 6. Annotate and write output ──────────────────────────────────────────────
sig_tiles$neuron_specific_overlap <- 
    seq_len(nrow(sig_tiles)) %in% tiles_with_neuron_signal

sig_tiles$composition_flag <- ifelse(
    sig_tiles$neuron_specific_overlap,
    "potentially_composition_driven",
    "composition_robust"
)

#write.table(sig_tiles, "sig_tiles_annotated.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(sig_tiles, "sig_tiles_annotated_q010.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# ── 7. Summary ────────────────────────────────────────────────────────────────
cat("\n── Compositional robustness summary ──\n")
cat("Composition-robust tiles:", sum(!sig_tiles$neuron_specific_overlap), "\n")
cat("  Hyper:", sum(!sig_tiles$neuron_specific_overlap & sig_tiles$meth.diff > 0), "\n")
cat("  Hypo: ", sum(!sig_tiles$neuron_specific_overlap & sig_tiles$meth.diff < 0), "\n")
cat("Potentially composition-driven tiles:", sum(sig_tiles$neuron_specific_overlap), "\n")
cat("  Hyper:", sum(sig_tiles$neuron_specific_overlap & sig_tiles$meth.diff > 0), "\n")
cat("  Hypo: ", sum(sig_tiles$neuron_specific_overlap & sig_tiles$meth.diff < 0), "\n")
cat("\nWritten to sig_tiles_annotated_q010.tsv\n")
