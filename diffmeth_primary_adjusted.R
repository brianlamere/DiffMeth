suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

# ── 1. Covariates from scMD deconvolution proportions ────────────────────────
prop <- fread("scMD_proportions.tsv")   # first col = sample (rowname), then 7 cell types
# scMD_proportions.tsv was written with row.names=TRUE, so first column is blank-headed
setnames(prop, 1, "sample")

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
prop <- prop[match(sample_order, prop$sample)]   # enforce order
stopifnot(all(prop$sample == sample_order))

neuron_frac <- prop$Inh + prop$Exc
oligo_lin   <- prop$Oligo + prop$OPC

covariates <- data.frame(
    neuron       = neuron_frac,
    oligo_lineage = oligo_lin
)
cat("Covariates (per sample):\n")
print(data.frame(sample=sample_order, cohort=c(rep("toxo+",6),rep("toxo-",3)),
                 neuron=round(neuron_frac,3), oligo_lineage=round(oligo_lin,3)))

# quick collinearity sanity check vs treatment
treatment <- c(0,0,0,0,0,0,1,1,1)   # 0=toxo+ (2024), 1=toxo- (2026)
cat("\nCorrelation of covariates with treatment (watch for ~|1|):\n")
cat("  neuron vs treatment:       ", round(cor(neuron_frac, treatment), 3), "\n")
cat("  oligo_lineage vs treatment:", round(cor(oligo_lin,  treatment), 3), "\n")

# ── 2. methylKit pipeline (identical to validated prior runs) ────────────────
sample_files <- as.list(paste0(sample_order, "_methylkit.txt"))
sample_ids   <- as.list(sample_order)

myobj <- methRead(location=sample_files, sample.id=sample_ids,
                  assembly="hg38", treatment=treatment,
                  context="CpG", mincov=10)
filtered     <- filterByCoverage(myobj, lo.count=10, hi.perc=99.9)
normalized   <- normalizeCoverage(filtered)
united_tiles <- tileMethylCounts(normalized) |> unite(min.per.group=3L)

# ── 3. Differential: unadjusted + primary adjusted (neuron + oligo-lineage) ──
dm_unadj <- calculateDiffMeth(united_tiles, overdispersion="MN", mc.cores=8)
dm_primary <- calculateDiffMeth(united_tiles, overdispersion="MN",
                                covariates=covariates, mc.cores=8)

# ── 4. Outputs ───────────────────────────────────────────────────────────────
# Full unadjusted dump (all ~86k tiles) — she filters herself
fwrite(getData(dm_unadj), "all_tiles_unadjusted.csv")

# Full primary-adjusted dump (all tiles)
fwrite(getData(dm_primary), "all_tiles_primary_adjusted.csv")

# Significant sets at q<0.05 for a quick comparison
sig_unadj   <- getMethylDiff(dm_unadj,   difference=10, qvalue=0.05)
sig_primary <- getMethylDiff(dm_primary, difference=10, qvalue=0.05)

key <- function(x) paste(getData(x)$chr, getData(x)$start, getData(x)$end)
ku <- key(sig_unadj); kp <- key(sig_primary)

cat("\n══ Significant tiles (q<0.05, |diff|>10%) ══\n")
cat("Unadjusted:         ", nrow(sig_unadj), "\n")
cat("Primary adjusted:   ", nrow(sig_primary), "\n")
cat("  survived adjustment:", sum(ku %in% kp), "\n")
cat("  lost after adjustment:", sum(!ku %in% kp), "\n")
cat("  new after adjustment:", sum(!kp %in% ku), "\n")

fwrite(getData(sig_unadj),   "sig_unadjusted.csv")
fwrite(getData(sig_primary), "sig_primary_adjusted.csv")

cat("\nWrote:\n")
cat("  all_tiles_unadjusted.csv        (all ~86k tiles, unadjusted)\n")
cat("  all_tiles_primary_adjusted.csv  (all tiles, neuron + oligo-lineage adjusted)\n")
cat("  sig_unadjusted.csv / sig_primary_adjusted.csv (q<0.05 subsets)\n")
