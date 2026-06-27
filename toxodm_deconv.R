#!/usr/bin/env Rscript
# ToxoDM cell-type deconvolution (phase 2, full 23-sample set)
# ---------------------------------------------------------------------------
# Method: NNLS (non-negative least squares) against the Lee single-cell WGBS
# reference signature (Lee et al. 2019, Nat Methods), as distributed in the scMD
# package. NNLS is the appropriate reference-based deconvolution method for
# methylation beta data; we use ONLY the scMD reference object, not its ensemble
# orchestration (which is non-reproducible on coordinate WGBS data and mixes in
# transcriptomics-oriented methods unsuited to methylation).
#
# Validation: the pipeline self-tests against the Guintivano sorted dataset
# (NeuN+/NeuN- known fractions) before running, asserting that NNLS recovers the
# neuron fraction (r > 0.99). This guards against a consistently-wrong stack
# (e.g. a mis-built BLAS) producing reproducible-but-inaccurate results -- a
# failure mode that re-run determinism alone cannot catch.
#
# Output: a 4-type proportion table (Astro, Micro, Neuron, Oligo) formed by
# summing OUTPUT fractions (Neuron = Inh + Exc). Summing output fractions is
# valid arithmetic (fractions are additive) and is the form validated against
# Guintivano. The unstable within-neuron Inh/Exc split is summed away. Endo and
# OPC are dropped (rare, unstable, not used as covariates). The 7-type table is
# also written for transparency/provenance.
#
# A settings manifest is exported alongside the proportions as a first-class
# deliverable, so any result is never an orphan from the parameters that made it.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(scMD)        # USED ONLY for the Lee WGBS reference object + Guintivano data
  library(nnls)
  library(data.table)
})

ALIGN  <- "/projects1/ToxoDM/working/align"
OUTDIR <- "/projects1/ToxoDM/reports"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# cell-type groupings for the 4-type collapse (sum of output fractions)
NEURON_TYPES <- c("Inh","Exc")        # summed -> Neuron
DROP_TYPES   <- c("Endo","OPC")       # rare/unstable, not covariates
# Astro, Micro, Oligo carried through as-is

# ---- helper: NNLS deconvolution of a bulk matrix against a reference -----------
# bulk: CpGs (rows, matched ids) x samples (cols);  ref: same CpG rows x celltypes
nnls_deconv <- function(bulk_mat, ref_mat) {
  common <- intersect(rownames(ref_mat), rownames(bulk_mat))
  refc  <- ref_mat[common, , drop = FALSE]
  bulkc <- bulk_mat[common, , drop = FALSE]
  out <- t(apply(bulkc, 2, function(col) {
    ok <- !is.na(col)
    fit <- nnls(refc[ok, , drop = FALSE], col[ok])
    fit$x / sum(fit$x)
  }))
  colnames(out) <- colnames(ref_mat)
  out
}

# =================================================================================
# STEP 0 -- VALIDATION GATE: NNLS must recover Guintivano known neuron fractions
# =================================================================================
# Guintivano is 450k (cg-probe ids) with NeuN+/NeuN- ground truth, so we validate
# against the 450k Lee reference. This validates the METHOD + STACK on the neuron
# axis (the high-stakes covariate). It does not test RRBS sparsity (Guintivano is
# full-coverage array) -- that remains a stated limitation, not tested here.
cat("== VALIDATION GATE: Guintivano neuron-fraction recovery ==\n")
data("Guintivano")        # Guintivano_bulk_sub, Guintivano_truefrac_sub
data("Lee_7ct_450850")    # Lee_sig_all (450k cg-id reference)

g_prop   <- nnls_deconv(Guintivano_bulk_sub, Lee_sig_all)
g_neuron <- rowSums(g_prop[, intersect(NEURON_TYPES, colnames(g_prop)), drop = FALSE])
g_truth  <- Guintivano_truefrac_sub[, "NeuN_pos"] /
            rowSums(Guintivano_truefrac_sub[, c("NeuN_pos","NeuN_neg"), drop = FALSE])
s <- intersect(names(g_neuron), names(g_truth))
g_r    <- cor(g_neuron[s], g_truth[s])
g_rmse <- sqrt(mean((g_neuron[s] - g_truth[s])^2))
cat(sprintf("  neuron-fraction recovery: r = %.4f, RMSE = %.4f (n=%d)\n",
            g_r, g_rmse, length(s)))

VALIDATION_R_MIN <- 0.99
if (is.na(g_r) || g_r < VALIDATION_R_MIN) {
  stop(sprintf(paste0("VALIDATION FAILED: neuron-fraction recovery r = %.4f < %.2f.\n",
       "  The deconvolution stack does not recover known fractions -- refusing to\n",
       "  produce sample proportions. Check R/BLAS/nnls/reference integrity."),
       g_r, VALIDATION_R_MIN))
}
cat(sprintf("  PASS (r >= %.2f). Proceeding to sample deconvolution.\n\n", VALIDATION_R_MIN))

# =================================================================================
# STEP 1 -- discover the 23 beds (follow symlinks for inherited CC)
# =================================================================================
beds <- list.files(ALIGN, pattern = "_CpG\\.bed$", full.names = TRUE, recursive = TRUE)
stopifnot(length(beds) == 23)

meta <- rbindlist(lapply(beds, function(p) {
  parts <- strsplit(sub(paste0("^", ALIGN, "/"), "", p), "/")[[1]]  # cohort/tissue/sample/file
  data.table(cohort = parts[1], tissue = parts[2], sample = parts[3], path = p)
}))
meta[, key := paste(cohort, tissue, sample, sep = "_")]
cat("Discovered", nrow(meta), "beds:\n")
print(meta[, .(cohort, tissue, sample)])

# =================================================================================
# STEP 2 -- build bulk beta matrix keyed by chr:pos
# =================================================================================
# bed cols: chr start end beta coverage  (0-based start = C position)
read_bed <- function(p) {
  b <- fread(p, col.names = c("chr","start","end","beta","coverage"))
  b[, coord := paste0(sub("^chr","",chr), ":", start)]   # match Lee_sig_all_WGBS rownames
  b[, .(coord, beta)]
}

cat("\nReading beds and assembling bulk matrix...\n")
bedlist <- lapply(meta$path, read_bed)
names(bedlist) <- meta$key

all_coords <- unique(unlist(lapply(bedlist, function(x) x$coord)))
bulk <- matrix(NA_real_, nrow = length(all_coords), ncol = nrow(meta),
               dimnames = list(all_coords, meta$key))
for (i in seq_len(nrow(meta))) {
  b <- bedlist[[i]]
  bulk[b$coord, meta$key[i]] <- b$beta
}
cat("Bulk matrix:", nrow(bulk), "CpGs x", ncol(bulk), "samples\n")

# =================================================================================
# STEP 3 -- restrict to the Lee WGBS signature CpGs, complete cases across samples
# =================================================================================
sig_coords <- rownames(Lee_sig_all_WGBS)
bulk_sig <- bulk[rownames(bulk) %in% sig_coords, , drop = FALSE]
bulk_cc  <- bulk_sig[rowSums(is.na(bulk_sig)) == 0, , drop = FALSE]
n_cc <- nrow(bulk_cc)
cat("Signature CpGs covered in ALL", ncol(bulk), "samples:", n_cc, "\n")

# =================================================================================
# STEP 4 -- NNLS deconvolution (7-type) against the Lee WGBS reference
# =================================================================================
cat("\nRunning NNLS deconvolution (Lee WGBS reference, 7-type)...\n")
ref7  <- Lee_sig_all_WGBS[rownames(bulk_cc), , drop = FALSE]   # n_cc x 7
prop7 <- nnls_deconv(bulk_cc, ref7)

# =================================================================================
# STEP 5 -- collapse to 4-type by SUMMING OUTPUT FRACTIONS (the validated form)
# =================================================================================
# Neuron = Inh + Exc (summed); Oligo, Astro, Micro as-is; Endo, OPC dropped.
# After dropping Endo/OPC the remaining 4 are renormalized to sum to 1, so the
# table is a clean 4-type composition the downstream tool reads with no further
# manipulation.
neuron <- rowSums(prop7[, intersect(NEURON_TYPES, colnames(prop7)), drop = FALSE])
prop4 <- data.frame(
  Astro_NNLS        = prop7[, "Astro"],
  Micro_NNLS        = prop7[, "Micro"],
  Neuron_total_NNLS = neuron,
  Oligo_NNLS        = prop7[, "Oligo"]
)
prop4 <- prop4 / rowSums(prop4)                 # renormalize after dropping Endo/OPC
prop4 <- round(prop4, 4)

# attach metadata, ordered to match meta
ord  <- match(rownames(prop4), meta$key)
out4 <- cbind(meta[ord, .(cohort, tissue, sample)], prop4)
out7 <- cbind(meta[match(rownames(prop7), meta$key), .(cohort, tissue, sample)],
              round(as.data.frame(prop7), 4))

cat("\n== 4-type proportions (Astro, Micro, Neuron, Oligo) ==\n")
print(out4, row.names = FALSE)

# =================================================================================
# STEP 6 -- write deliverables: 4-type table, 7-type table, settings manifest
# =================================================================================
f4 <- file.path(OUTDIR, "nnls_proportions_4type.tsv")
f7 <- file.path(OUTDIR, "nnls_proportions_7type.tsv")
fm <- file.path(OUTDIR, "nnls_proportions_SETTINGS.txt")
fwrite(out4, f4, sep = "\t")
fwrite(out7, f7, sep = "\t")

# settings manifest: self-describing, 2077-legible, shipped WITH the results
manifest <- c(
  "ToxoDM cell-type deconvolution -- settings manifest",
  paste0("generated:            ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "METHOD",
  "  algorithm:            NNLS (non-negative least squares), normalized to sum 1",
  "  reference:            Lee_sig_all_WGBS (Lee et al. 2019 Nat Methods, via scMD pkg)",
  "  reference n_celltypes:7 (Astro Micro Endo Oligo OPC Inh Exc)",
  paste0("  scMD pkg version:     ", as.character(packageVersion("scMD"))),
  paste0("  nnls pkg version:     ", as.character(packageVersion("nnls"))),
  paste0("  R version:            ", R.version.string),
  "",
  "INPUT",
  paste0("  bulk source:          ", ALIGN, "/<cohort>/<tissue>/<sample>/*_CpG.bed"),
  paste0("  n_samples:            ", nrow(meta)),
  paste0("  signature CpGs (cc):  ", n_cc, " (covered in all samples)"),
  "  bulk values:          biscuit CpG beta (strand-collapsed), no coverage cutoff here",
  "",
  "4-TYPE COLLAPSE (by summing OUTPUT fractions; fractions are additive)",
  paste0("  Neuron_total_NNLS = ", paste(NEURON_TYPES, collapse = " + ")),
  "  Astro_NNLS, Micro_NNLS, Oligo_NNLS = corresponding 7-type fractions as-is",
  paste0("  dropped:              ", paste(DROP_TYPES, collapse = ", "),
         " (rare/unstable, not used as covariates)"),
  "  renormalized:         remaining 4 types scaled to sum to 1 per sample",
  "",
  "VALIDATION (stack + method sanity gate, run before sample deconvolution)",
  "  set:                  Guintivano sorted (NeuN+/NeuN- known fractions, 450k)",
  "  axis validated:       neuron fraction (Neuron = Inh+Exc summed)",
  paste0("  result:               r = ", sprintf("%.4f", g_r),
         ", RMSE = ", sprintf("%.4f", g_rmse), " (n=", length(s), ")"),
  paste0("  gate threshold:       r >= ", VALIDATION_R_MIN, " (PASS)"),
  "  NOT validated:        within-glia split (Astro/Micro/Oligo), RRBS sparsity",
  "                        (Guintivano is full-coverage array, not RRBS)",
  "",
  "NOTES",
  "  - scMD ensemble orchestration intentionally NOT used: non-reproducible on",
  "    coordinate WGBS data (order/state-dependent), and mixes transcriptomics",
  "    methods unsuited to methylation. Only the validated Lee reference is used.",
  "  - Absolute cross-tissue proportions carry method bias; values are intended",
  "    as within-tissue relative covariates (bias ~constant within tissue cancels)."
)
writeLines(manifest, fm)

cat("\nWrote:\n  ", f4, "\n  ", f7, "\n  ", fm, "\n")
