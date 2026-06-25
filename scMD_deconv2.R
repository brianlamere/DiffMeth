suppressPackageStartupMessages(library(scMD))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
beds <- lapply(sample_order, function(s) {
    b <- fread(file.path(s, paste0(s,"_GRCh38_CpG.bed")),
               col.names=c("chr","start","end","beta","coverage"))
    b[, sig_coord := paste0(sub("^chr","",chr), ":", start)]
    b[, .(sig_coord, beta)]
})
names(beds) <- sample_order
all_coords <- unique(unlist(lapply(beds, function(x) x$sig_coord)))
bulk <- matrix(NA_real_, nrow=length(all_coords), ncol=length(sample_order),
               dimnames=list(all_coords, sample_order))
for (s in sample_order) {
    m <- match(beds[[s]]$sig_coord, all_coords)
    bulk[m, s] <- beds[[s]]$beta
}

# Restrict to signature CpGs covered in all 9 samples (complete cases)
sig_coords <- rownames(Lee_sig_all_WGBS)
bulk_sig  <- bulk[rownames(bulk) %in% sig_coords, , drop=FALSE]
bulk_cc   <- bulk_sig[rowSums(is.na(bulk_sig))==0, , drop=FALSE]
cat("Complete-case signature CpGs feeding deconvolution:", nrow(bulk_cc), "\n")

# Run scMD WITHOUT Houseman (which forces 450k probe names via minfi).
# NNLS and RPC are coordinate-agnostic matrix methods.
result <- tryCatch(
    scMD(bulk = bulk_cc, bulk_type = "WGBS",
         dmet_list = c("NNLS","RPC","CIBERSORT","EPIC","FARDEEP","DCQ","ICeDT")),
    error = function(e) { cat("scMD error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(result)) {
    cat("\n══ scMD estimated proportions (ensemble) ══\n")
    print(round(result$scMD_p, 4))
    write.table(round(result$scMD_p,4), "scMD_proportions.tsv", sep="\t", quote=FALSE)
}
