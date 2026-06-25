suppressPackageStartupMessages(library(scMD))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

# scMD() wants a bulk methylation matrix: rows = CpGs, cols = samples,
# values = beta (methylation fraction). It matches against its signature
# internally. For "450k_or_850k" it expects cg-probe rownames; for "WGBS"
# it should accept coordinate rownames matching its WGBS signature.

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")

# Build a bulk beta matrix keyed by chr:pos (to match Lee_sig_all_WGBS coords)
beds <- lapply(sample_order, function(s) {
    b <- fread(file.path(s, paste0(s,"_GRCh38_CpG.bed")),
               col.names=c("chr","start","end","beta","coverage"))
    # signature coords are like "5:34189635" -> strip chr, use start
    b[, sig_coord := paste0(sub("^chr","",chr), ":", start)]
    b[, .(sig_coord, beta)]
})
names(beds) <- sample_order

# Union of all covered coords, then assemble matrix (NA where a sample lacks coverage)
all_coords <- unique(unlist(lapply(beds, function(x) x$sig_coord)))
cat("Total unique covered CpGs across samples:", length(all_coords), "\n")

bulk <- matrix(NA_real_, nrow=length(all_coords), ncol=length(sample_order),
               dimnames=list(all_coords, sample_order))
for (s in sample_order) {
    m <- match(beds[[s]]$sig_coord, all_coords)
    bulk[m, s] <- beds[[s]]$beta
}

# How much overlaps the scMD WGBS signature?
sig_coords <- rownames(Lee_sig_all_WGBS)
covered_in_sig <- sum(rownames(bulk) %in% sig_coords)
cat("Of those, in scMD WGBS signature:", covered_in_sig, "\n")
cat("Signature size:", length(sig_coords), "\n")
cat("Fraction of signature covered:", round(covered_in_sig/length(sig_coords), 4), "\n")

# scMD with complete-case requirement will be brutal on sparse RRBS;
# restrict bulk to signature CpGs and see how many are covered in ALL samples
bulk_sig <- bulk[rownames(bulk) %in% sig_coords, , drop=FALSE]
complete <- bulk_sig[rowSums(is.na(bulk_sig))==0, , drop=FALSE]
cat("Signature CpGs covered in ALL 9 samples:", nrow(complete), "\n\n")

# Attempt the deconvolution
cat("Attempting scMD deconvolution (WGBS mode)...\n")
result <- tryCatch(
    scMD(bulk = bulk, bulk_type = "WGBS"),
    error = function(e) { cat("scMD error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(result)) {
    cat("\n══ scMD estimated proportions ══\n")
    print(round(result$scMD_p, 4))
    write.table(round(result$scMD_p,4), "scMD_proportions.tsv",
                sep="\t", quote=FALSE)
}
