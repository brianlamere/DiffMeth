#!/usr/bin/env Rscript
# Convert all 23 ToxoDM biscuit CpG beds -> methylKit input format.
# Walks working/align/<cohort>/<tissue>/<sample>/*_CpG.bed (follows symlinks for
# inherited CC). Writes <sample-key>_methylkit.txt into working/methylkit/.
# Pure format reshuffle (cheap); regenerated fresh so nothing depends on phase-1 files.

ALIGN  <- "/projects1/ToxoDM/working/align"
OUTDIR <- "/projects1/ToxoDM/working/methylkit"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

convert_biscuit_to_methylkit <- function(bed_file, out_file) {
  df <- read.table(bed_file, header = FALSE, sep = "\t",
                   col.names = c("chr","start","end","beta","coverage"))
  df <- df[df$coverage > 0, ]                       # drop zero-coverage (defensive)
  out <- data.frame(
    chrBase  = paste(df$chr, df$start, sep = "."),
    chr      = df$chr,
    base     = df$start,
    strand   = "F",                                 # biscuit beds are strand-collapsed CpG
    coverage = df$coverage,
    freqC    = round(df$beta * 100, 4),
    freqT    = round((1 - df$beta) * 100, 4)
  )
  write.table(out, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

beds <- list.files(ALIGN, pattern = "_CpG\\.bed$", full.names = TRUE, recursive = TRUE)
stopifnot(length(beds) == 23)

for (p in beds) {
  parts  <- strsplit(sub(paste0("^", ALIGN, "/"), "", p), "/")[[1]]  # cohort/tissue/sample/file
  key    <- paste(parts[1], parts[2], parts[3], sep = "_")
  outf   <- file.path(OUTDIR, paste0(key, "_methylkit.txt"))
  message("Converting ", key)
  convert_biscuit_to_methylkit(p, outf)
}
cat("\nWrote", length(beds), "methylKit files to", OUTDIR, "\n")
