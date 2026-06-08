library(methylKit)

convert_biscuit_to_methylkit <- function(bed_file, out_file) {
  df <- read.table(bed_file, header=FALSE, sep="\t",
                   col.names=c("chr","start","end","beta","coverage"))
  
  # Drop zero coverage rows (shouldn't be any, but just in case)
  df <- df[df$coverage > 0, ]
  
  # Build methylKit format
  out <- data.frame(
    chrBase  = paste(df$chr, df$start, sep="."),
    chr      = df$chr,
    base     = df$start,
    strand   = "F",
    coverage = df$coverage,
    freqC    = round(df$beta * 100, 4),
    freqT    = round((1 - df$beta) * 100, 4)
  )
  
  write.table(out, file=out_file, sep="\t", quote=FALSE, row.names=FALSE)
}

# Run for all 9 samples
samples <- list(
  list(bed="CA192/CA192_GRCh38_mq30_CpG.bed", out="CA192_methylkit_mq30.txt"),
  list(bed="CA346/CA346_GRCh38_mq30_CpG.bed", out="CA346_methylkit_mq30.txt"),
  list(bed="CB239/CB239_GRCh38_mq30_CpG.bed", out="CB239_methylkit_mq30.txt"),
  list(bed="CC249/CC249_GRCh38_mq30_CpG.bed", out="CC249_methylkit_mq30.txt"),
  list(bed="CE167/CE167_GRCh38_mq30_CpG.bed", out="CE167_methylkit_mq30.txt"),
  list(bed="CE234/CE234_GRCh38_mq30_CpG.bed", out="CE234_methylkit_mq30.txt"),
  list(bed="LG30/LG30_GRCh38_mq30_CpG.bed",   out="LG30_methylkit_mq30.txt"),
  list(bed="LG31/LG31_GRCh38_mq30_CpG.bed",   out="LG31_methylkit_mq30.txt"),
  list(bed="LG52/LG52_GRCh38_mq30_CpG.bed",   out="LG52_methylkit_mq30.txt")
)

treatment = c(0,0,0,0,0,0,1,1,1)

for (s in samples) {
  message("Converting ", s$bed)
  convert_biscuit_to_methylkit(s$bed, s$out)
}
