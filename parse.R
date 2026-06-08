library(methylKit)

# Build file list from your BED outputs
# methylKit can read bismark coverage format natively;
# for BISCUIT beds you'll want a small conversion step

sample_files <- list(
  "CA192_methylkit.txt", "CA346_methylkit.txt", "CB239_methylkit.txt", "CC249_methylkit.txt",
  "CE167_methylkit.txt", "CE234_methylkit.txt", "LG30_methylkit.txt", "LG31_methylkit.txt", "LG52_methylkit.txt"
)

myobj <- methRead(
  location = sample_files,
  sample.id = list("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52"),
  assembly = "hg38",
  treatment = c(0,0,0,0,0,0,1,1,1),  # 2024=0, 2026=1
  context = "CpG",
  mincov = 10  # adjust; 10 is reasonable for RRBS
)

# Filter for coverage
filtered <- filterByCoverage(myobj, lo.count=10, hi.perc=99.9)

# Normalize (tiles reads across samples)
normalized <- normalizeCoverage(filtered)

# Unite — only CpGs covered in all samples, or adjust min.per.group
united <- unite(normalized, min.per.group=3L)

# Differential methylation
dm_cpg <- calculateDiffMeth(united,
  overdispersion="MN",  # better for biological replicates
  test="Chisq",
  mc.cores=8
)

dm_tiles <- tileMethylCounts(normalized) |>
  unite(min.per.group=3L) |>
  calculateDiffMeth(overdispersion="MN", mc.cores=8)

# Get significant differentially methylated tiles
sig_tiles <- getMethylDiff(dm_tiles, difference=10, qvalue=0.05)
nrow(sig_tiles)  # how many pass threshold

# Or look at the full distribution
hist(dm_tiles$qvalue, breaks=50)
hist(dm_tiles$meth.diff, breaks=100)

# How many tiles go each direction among significant hits
sig_tiles_hyper <- getMethylDiff(dm_tiles, difference=10, qvalue=0.05, type="hyper")  # higher in 2026
sig_tiles_hypo  <- getMethylDiff(dm_tiles, difference=10, qvalue=0.05, type="hypo")   # lower in 2026

nrow(sig_tiles_hyper)
nrow(sig_tiles_hypo)

# Summary stats to quantify the asymmetry
summary(dm_tiles$meth.diff)

# How many tiles exceed various thresholds
table(cut(dm_tiles$meth.diff, breaks=c(-Inf,-25,-10,10,25,Inf)))

# Check if the asymmetry holds among your most significant hits only
sig <- getMethylDiff(dm_tiles, difference=10, qvalue=0.05)
nrow(sig)
table(sig$meth.diff > 0)  # hyper vs hypo among significant

sig_coords <- getData(sig)[, c("chr","start","end")]

# Coverage is stored in the united object
# percMethylation gives you the matrix, but for coverage:
cov_matrix <- getCoverageMatrix(united)  # if this method exists in your version

# More reliably:
coverage_data <- getData(united)
# This gives you all the raw columns including coverage per sample
head(coverage_data)
