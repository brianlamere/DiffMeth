suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

setwd("/projects2/DiffMeth")

# ══ 1. Build per-sample oligodendrocyte covariate values ════════════════════
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

markers <- c("cg16515546", "cg09837648")
mc <- get_hg38_coords(markers, anno_hg19, anno_v2, anno_v2_base, chain)

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234",
                  "LG30","LG31","LG52")
sample_beds <- setNames(
    file.path(sample_order, paste0(sample_order, "_GRCh38_CpG.bed")),
    sample_order)

# Pull beta + depth for each marker in each sample
#get_marker_vals <- function(coord) {
#    sapply(sample_order, function(s) {
#        bed <- fread(sample_beds[[s]],
#                     col.names=c("chr","start","end","beta","coverage"))
#        bed$coord <- paste(bed$chr, bed$start, sep=":")
#        row <- bed[coord == ..coord]   # note: ..coord to reference outer var
#        if (nrow(row) > 0) c(beta=row$beta[1], cov=row$coverage[1])
#        else c(beta=NA, cov=NA)
#    })
#}

# (simpler explicit loop to avoid data.table scoping subtleties)
get_marker_vals <- function(target_coord) {
    out <- matrix(NA, nrow=2, ncol=length(sample_order),
                  dimnames=list(c("beta","cov"), sample_order))
    for (s in sample_order) {
        bed <- fread(sample_beds[[s]],
                     col.names=c("chr","start","end","beta","coverage"))
        bed$coord <- paste(bed$chr, bed$start, sep=":")
        row <- bed[bed$coord == target_coord, ]
        if (nrow(row) > 0) {
            out["beta", s] <- row$beta[1]
            out["cov",  s] <- row$coverage[1]
        }
    }
    out
}

m1 <- get_marker_vals(mc$coord[mc$probe_id=="cg16515546"])
m2 <- get_marker_vals(mc$coord[mc$probe_id=="cg09837648"])

# Covariate A: cg16515546 alone
covA <- data.frame(oligo = as.numeric(m1["beta", ]))

# Covariate B: coverage-weighted composite of both markers
composite <- sapply(sample_order, function(s) {
    b1 <- m1["beta", s]; c1 <- m1["cov", s]
    b2 <- m2["beta", s]; c2 <- m2["cov", s]
    vals <- c(b1, b2); ws <- c(c1, c2)
    ok <- !is.na(vals) & !is.na(ws)
    if (any(ok)) sum(vals[ok]*ws[ok]) / sum(ws[ok]) else NA
})
covB <- data.frame(oligo = as.numeric(composite))

cat("Covariate A (cg16515546 alone):\n"); print(round(covA$oligo, 3))
cat("\nCovariate B (composite):\n");      print(round(covB$oligo, 3))

# ══ 2. methylKit pipeline ════════════════════════════════════════════════════
sample_files <- as.list(file.path(paste0(sample_order, "_methylkit.txt")))
sample_ids   <- as.list(sample_order)
treatment    <- c(0,0,0,0,0,0,1,1,1)

myobj <- methRead(location=sample_files, sample.id=sample_ids,
                  assembly="hg38", treatment=treatment,
                  context="CpG", mincov=10)
filtered   <- filterByCoverage(myobj, lo.count=10, hi.perc=99.9)
normalized <- normalizeCoverage(filtered)
united_tiles <- tileMethylCounts(normalized) |> unite(min.per.group=3L)

# ══ 3. Three differential runs ═══════════════════════════════════════════════
dm_unadj <- calculateDiffMeth(united_tiles, overdispersion="MN", mc.cores=8)
dm_covA  <- calculateDiffMeth(united_tiles, overdispersion="MN",
                              covariates=covA, mc.cores=8)
dm_covB  <- calculateDiffMeth(united_tiles, overdispersion="MN",
                              covariates=covB, mc.cores=8)

sig_unadj <- getMethylDiff(dm_unadj, difference=10, qvalue=0.05)
sig_covA  <- getMethylDiff(dm_covA,  difference=10, qvalue=0.05)
sig_covB  <- getMethylDiff(dm_covB,  difference=10, qvalue=0.05)

#adding this looser threshold for exploration/validation, to be reviewed "manually"
sig_unadjq1 <- getMethylDiff(dm_unadj, difference=10, qvalue=0.1)
sig_covAq1  <- getMethylDiff(dm_covA,  difference=10, qvalue=0.1)
sig_covBq1  <- getMethylDiff(dm_covB,  difference=10, qvalue=0.1)

write.table(getData(sig_unadjq1), "sig_unadj_q10.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(getData(sig_covAq1),  "sig_adj_cg16515546_q10.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(getData(sig_covBq1),  "sig_adj_composite_q10.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# ══ 4. Comparison ════════════════════════════════════════════════════════════
key <- function(x) paste(getData(x)$chr, getData(x)$start, getData(x)$end)
ku <- key(sig_unadj); kA <- key(sig_covA); kB <- key(sig_covB)

cat("\n══ Significant tile counts ══\n")
cat("Unadjusted:                ", nrow(sig_unadj), "\n")
cat("Adjusted (cg16515546):     ", nrow(sig_covA), "\n")
cat("Adjusted (composite):      ", nrow(sig_covB), "\n")

cat("\n══ vs unadjusted baseline ══\n")
cat("cg16515546 - survived:", sum(ku %in% kA),
    " lost:", sum(!ku %in% kA), " new:", sum(!kA %in% ku), "\n")
cat("composite  - survived:", sum(ku %in% kB),
    " lost:", sum(!ku %in% kB), " new:", sum(!kB %in% ku), "\n")

cat("\n══ covA vs covB agreement ══\n")
cat("In both adjusted sets:", sum(kA %in% kB), "\n")
cat("covA only:", sum(!kA %in% kB), "  covB only:", sum(!kB %in% kA), "\n")

# ══ 5. Save ══════════════════════════════════════════════════════════════════
write.table(getData(sig_unadj), "sig_unadj.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(getData(sig_covA),  "sig_adj_cg16515546.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(getData(sig_covB),  "sig_adj_composite.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\nWritten: sig_unadj.tsv, sig_adj_cg16515546.tsv, sig_adj_composite.tsv\n")
