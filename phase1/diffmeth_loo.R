suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(data.table))
setwd("/projects2/DiffMeth")

prop <- fread("scMD_proportions.tsv"); setnames(prop, 1, "sample")
sample_order <- c("CA192","CA346","CB239","CC249","CE167","CE234","LG30","LG31","LG52")
prop <- prop[match(sample_order, prop$sample)]
treatment_full <- c(0,0,0,0,0,0,1,1,1)
oligo_lin_full <- prop$Oligo + prop$OPC

sample_files_full <- paste0(sample_order, "_methylkit.txt")

# Run one differential on a given subset of sample indices.
# adjust = NULL (unadjusted) or "oligolin"
run_one <- function(keep_idx, adjust) {
    ids   <- as.list(sample_order[keep_idx])
    files <- as.list(sample_files_full[keep_idx])
    trt   <- treatment_full[keep_idx]
    # need >=2 per group to compare; min.per.group set accordingly
    n_ctrl <- sum(trt==1); n_case <- sum(trt==0)
    mpg <- min(3L, n_ctrl, n_case)            # drop to 2 only when forced
    obj <- methRead(location=files, sample.id=ids, assembly="hg38",
                    treatment=trt, context="CpG", mincov=10)
    f <- filterByCoverage(obj, lo.count=10, hi.perc=99.9)
    n <- normalizeCoverage(f)
    ut <- tileMethylCounts(n) |> unite(min.per.group=as.integer(mpg))
    cov <- if (is.null(adjust)) NULL else
           data.frame(oligo_lineage = oligo_lin_full[keep_idx])
    dm <- calculateDiffMeth(ut, overdispersion="MN",
                            covariates=cov, mc.cores=8)
    #list(sig = getMethylDiff(dm, difference=10, qvalue=0.05), mpg = mpg)
    list(sig = getMethylDiff(dm, difference=10, qvalue=0.1), mpg = mpg)
}

key <- function(x) paste(getData(x)$chr, getData(x)$start, getData(x)$end)

loo <- function(adjust, label) {
    cat("\n══ LOO:", label, "══\n")
    # baseline (all 9)
    base <- run_one(1:9, adjust)
    base_key <- key(base$sig)
    cat(sprintf("Baseline (all 9, min.per.group=%d): %d sig tiles\n\n",
                base$mpg, length(base_key)))
    cat(sprintf("%-8s %4s %6s %8s %6s %5s\n",
                "dropped","mpg","n_sig","retained","lost","new"))
    cat(strrep("-",46),"\n")
    res <- list()
    for (i in 1:9) {
        r <- run_one(setdiff(1:9, i), adjust)
        k <- key(r$sig)
        cat(sprintf("%-8s %4d %6d %8d %6d %5d\n",
            sample_order[i], r$mpg, length(k),
            sum(base_key %in% k), sum(!base_key %in% k), sum(!k %in% base_key)))
        res[[sample_order[i]]] <- k
    }
    list(base=base_key, drops=res)
}

unadj_loo <- loo(NULL,      "Unadjusted")
oligo_loo <- loo("oligolin","Oligo-lineage adjusted")

# Save which baseline tiles are lost when each sample is dropped
save_loo <- function(loo_res, tag) {
    bk <- loo_res$base
    dt <- data.table(tile = bk)
    for (s in names(loo_res$drops)) dt[[paste0("drop_",s)]] <- bk %in% loo_res$drops[[s]]
    dt[, n_drops_retaining := rowSums(.SD), .SDcols=patterns("^drop_")]
    fwrite(dt, paste0("loo_", tag, "q0.1.csv"))
}
save_loo(unadj_loo, "unadjusted")
save_loo(oligo_loo, "oligolin")
cat("\nWrote loo_unadjusted.csv and loo_oligolin.csv\n")
