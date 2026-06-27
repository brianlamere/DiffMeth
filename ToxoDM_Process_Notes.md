# ToxoDM / RRBS DiffMeth — Process Notes (technical struggles & decisions)

Notes-level record of the technical problems hit and the reasoning behind the
choices, across all phases of the toxo RRBS work (phase-1 9-sample cortex →
phase-2 23-sample multi-tissue). Not comprehensive — a skeleton to mod with the
activities I didn't narrate. Roughly chronological within each area.

---

## Alignment / methylation extraction

**Why biscuit (over bismark / bwa-meth).** Prior analyst's data had aligned
poorly to GRCh38. Chose biscuit on a bwa-mem2 core: bisulfite-aware, handles
soft-clipping well (matters for RRBS adapter remnants), recovers more reads on
difficult data. bismark's two-genome approach has higher memory overhead and is
stricter on mismatches; bwa-meth is fast but less maintained and has known
issues with some RRBS read structures. Used GRCh38_no_alt_analysis_set — the
no-alt set avoids multi-mapping inflation from alt contigs.

**No deduplication (the key RRBS call).** RRBS fragments are enzymatically
defined at MspI cut sites, so reads legitimately pile up at identical
start/end coordinates. With no UMIs, you cannot distinguish PCR duplicates from
biological stacking — deduplicating systematically wipes coverage, worst at the
high-methylation regions that are most interesting. A prior consultant ran
aggressive dedup on the 2024 data and zeroed it to ~no coverage (the data Sarah
couldn't use). First pass here did NOT deduplicate at all; that alone produced
usable data. (Standard-correct practice; don't owe a justification for it.)

**MQ filtering is already handled by biscuit.** Ran a parallel arm with
`samtools view -q 30` pre-filtering to test whether low-MQ reads drove results.
Output was byte-identical to the unfiltered arm. Reason: biscuit pileup applies
its own internal minimum mapping quality, default **MQ40** — stricter than the
samtools MQ30 pre-filter, so both fed identical reads into methylation calling.
Upshot: results were high-confidence (MQ40+) all along; the MQ30 arm added no
information and was dropped. (A looser `-m 30`/`-m 20` arm would only add noise
biscuit already judged ambiguous — not worth it except for a targeted
low-coverage locus near a repeat/segdup.)

**Pipeline shape.** trim_galore --rrbs --paired → biscuit align
GRCh38_no_alt → samtools sort → biscuit pileup (default MQ40) → vcf2bed -t cg →
convert to methylKit txt. methylKit over RNBeads for the control it gives on a
small asymmetric design with potential batch effects.

**RRBS QC normal characteristics** (so they don't look like errors): ~44.7%
base loss from Trim Galore is expected for short MspI fragments; duplicate
over-flagging is expected from deterministic MspI cutting; MAPQ distribution is
trimodal (multimappers / low-confidence / MQ60 high-confidence). CA346 is a
known thin/low-library donor — consistently lower depth, not a pipeline fault.
2026 LG samples run at higher depth than 2024 (cohort difference).

---

## Gentoo papercut: trim_galore vs python-exec2c

After a system update, trim_galore 0.6.11 died with "No Python detected"
despite cutadapt working fine. Cause: trim_galore reads cutadapt's shebang
(`#!/usr/bin/python-exec2c` on Gentoo), runs `python-exec2c --version`, gets
back `python-exec 2.4.10` which matches neither `Python 3.*` nor `Python 2.*`,
and bails. cutadapt itself was fine (built against python3.13, the active slot).
Fix: a portage post-install hook (`/etc/portage/bashrc`,
`post_pkg_postinst`) that rewrites the cutadapt shebang to `#!/usr/bin/python3`
after each rebuild (`/usr/bin/python3` follows the eselect-managed slot, so it's
forward-compatible). NOTE: trim_galore v2.x (Rust rewrite, GA May 2026) drops
the Python/cutadapt runtime dependency entirely — single static binary — which
makes this papercut disappear permanently if/when we migrate.

---

## Deconvolution — the long arc (HiBED → scMD → NNLS)

The composition-confound question: bulk RRBS methylation is a weighted average
over cell types. If composition differs between groups, apparent differential
methylation can be a composition shift, not biology. So we need per-sample
cell-type proportions to use as covariates.

**HiBED — failed: sparse-not-hitting-sparse.** HiBED uses array (450k) marker
probes. RRBS only covers ~5–10% of CpGs (those near MspI sites). The
informative deconvolution markers and the RRBS-covered sites are both sparse,
and they didn't overlap enough — sparse marker set tied to sparse coverage.
Only oligodendrocyte signal came through at the cohort-pool level; per-sample
resolution wasn't recoverable. Abandoned.

**scMD — extended debugging, ultimately abandoned the orchestration.** scMD is
a wrapper that runs EnsDeconv → multiple third-party deconvolution tools and
averages them ("ensemble"). It ships a Lee et al. 2019 WGBS reference
(coordinate-keyed, 7 brain cell types) that DID match our coordinate bulk
(~3,747 CpGs covered in all 23 samples). The reference is sound; the
orchestration is not. Characterized failure modes:

- **Non-reproducible / state-dependent.** The "good" proportions table (the one
  that would have ended the project) could only be produced from an accumulated
  warm R session, never cold from Rscript. Traced it: omitting `bulk_type`
  ran the default (450k) path first, loading state, then re-running with
  `bulk_type="WGBS"` drew on accumulated `phat_all` containing BOTH beta (valid)
  and Mval (all-NaN on RRBS, which is full of beta=0/1) entries plus duplicated
  NNLS. The table was a composite of multiple execution paths' leftovers, not a
  clean computation. A cold run either crashes at assembly or gives different
  numbers depending on path/order. The overlap proved the numbers came from the
  WGBS reference regardless of the bulk_type passed (450k overlap = 0, WGBS = 3747).

- **Assembly crash on partial failure.** `x[["a"]][["p_hat"]][[1]]` subscript-
  out-of-bounds when any method fails and leaves an empty p_hat. The ensemble
  can't survive its own partial failures — crashes even on scMD's OWN Guintivano
  example data, cold.

- **Houseman = hard halt.** Routes through minfi
  makeGenomicRatioSetFromMatrix, which hard-requires 450k cg-probe rownames —
  instant halt on coordinate data, before any save. Must be excluded.

- **RNA tools on methylation.** Most ensemble methods (CIBERSORT, EPIC,
  FARDEEP, DCQ, ICeDT) are transcriptomics deconvolution tools, built for gene
  expression / immune cells, misapplied to methylation. Produced incoherent
  per-sample estimates (one BG sample 47% neuron, the next 0%; oligo ranging
  8–64% across same-tissue samples). Only NNLS (and Houseman, which can't run)
  are methylation-appropriate. The "ensemble" dilutes the one good method (NNLS)
  with garbage from five wrong ones.

- **nmrk default bug (filed: github.com/randel/scMD/issues/2).** Both scMD()
  and sc_MD_deconv() declare `nmrk = 100` in their signatures and docs, yet the
  default execution runs at 50 (visible in the n_markers field of saved per-
  method records). Explicit values ARE honored (tested 50/100/500/5000; field
  updates, proportions change). So the documented/coded default of 100 doesn't
  reach the marker-selection code; some downstream default of 50 wins — and the
  real value lives in an unexported internal, invisible even in the source.
  Marker count is genuinely sensitive (per De Ridder 2024 Nat Comms benchmark,
  s41467-024-48466-z): changing it shifts the Inh/Exc split — more evidence the
  subtype split is a knob-dependent artifact, not a stable measurement. (Higher
  nmrk pushes the split toward NNLS's all-Inh corner; the plausible-looking Exc
  values in the warm-session table were partly a low-marker artifact.)

- **get_sig default ct_ind is 4 types** (Astro/Micro/Neuro/Oligo) but the
  shipped WGBS .rda has no native "Neuro" column — confirms the authors'
  intended resolution is 4-type and subtypes were WIP. The reference's own
  construction/aggregation code is opaque (Lee primary data is reputable, marker
  selection is transparent, but the aggregation foundry is unseen).

Decision: keep ONLY the validated Lee reference object; abandon all scMD
orchestration. Use NNLS directly.

**NNLS — the keeper.** NNLS (non-negative least squares) is the appropriate
reference-based method for methylation beta data. It's deterministic, runs cold
from Rscript, and there's no argument chain / inner default / session state to
contaminate it. Behavior: sparse solutions — it zeros out collinear pairs
(Exc→0, dumping neuronal signal into Inh; OPC→0 into Oligo). So the SUBTYPE
splits (Inh vs Exc, OPC vs Oligo) are not reliably estimable from sparse RRBS,
but the AGGREGATES are robust: neuron-total (Inh+Exc), oligo-lineage,
Micro, Astro. Collapse to 4 types by SUMMING OUTPUT FRACTIONS (valid arithmetic;
avoids the reference-collapse problem of averaging reference columns).

---

## Validation (the part that makes it defensible)

**Guintivano ground-truth check.** Guintivano = 12 sorted brain samples
(6 donors × NeuN+/NeuN−), 450k array, neuron fraction known by cell sorting.
Ran our NNLS process on Guintivano's bulk, summed Inh+Exc, scored vs NeuN_pos:
**r = 0.9991, RMSE = 0.168.** Near-perfect discrimination of neuron vs glia.
This validates the PROCESS (and the whole R/BLAS/nnls/reference stack) on
samples with a known answer, then we apply the same process to our unknowns.
Key framings:
- Validates neuron/non-neuron axis ONLY (all Guintivano sorted on). Does NOT
  validate the within-glia split or RRBS sparsity (Guintivano is full-coverage
  array). Stated as limitations, not claimed.
- The systematic compression bias (pure neurons read ~0.88 not 1.0; pure glia
  ~0.20 not 0) is a near-CONSTANT offset → cancels in within-tissue relative
  use (which is all the covariate does). Absolute cross-tissue proportions are
  NOT trustworthy from either method; relative within-tissue use is.
- The scMD ensemble CRASHED even on Guintivano — couldn't be scored at all.

**Aggregate-vs-known-tissue-biology as a second validation axis.** Composition
can't QC individual samples (autopsy dissection of e.g. basal ganglia inherently
varies white/grey — "where the ice cream scoop landed"). But aggregate over many
same-tissue samples should converge toward the tissue's true composition (law of
large numbers). NNLS aggregate: cortex ~30–36% neuron (plausible), BG ~29% (too
high — true BG is ~5–10%). Warm-scMD aggregate: BG ~6.5% (plausible), cortex
~11% (too low). So the two methods carry OPPOSITE absolute biases. NNLS is right
for cortex, scMD-ish for BG — but BG is never used as a neuron covariate, and
within-tissue the bias cancels, so NNLS stands. (LG30/LG52 at ~41% neuron =
plausibly clean grey-matter dissection, top-of-range but not alarming; don't
over-read the exact value — a few points could be compression bias.)

**USMC range analogy (why consistent-but-biased beats inconsistent-but-
sometimes-accurate).** Teach the non-hunters to group tight first, THEN adjust
the sights — a consistent-but-offset shooter is correctable; an inconsistent-
but-sometimes-accurate one is a lost cause, no systematic offset to dial out.
NNLS = tight group, sights adjusted via Guintivano. scMD ensemble = scattered,
can't be sighted in. Re-run determinism catches non-reproducibility; Guintivano
catches a consistently-WRONG stack (mis-built BLAS) that re-run determinism
can't. Built the Guintivano check in as a GATE: the deconv script refuses to
run (stop, r < 0.99) if the stack stops recovering known truth.

---

## Final pipeline shape (phase 2)

- Deconv: NNLS vs Lee WGBS reference, 7-type → collapse to 4 (Astro, Micro,
  Neuron=Inh+Exc, Oligo) by summing output fractions, renormalized. Guintivano
  self-test gates the run. Writes 4type + 7type tables + SETTINGS manifest.
- Differential battery (Sarah's plan): per comparison, four models —
  model_0 (~status), model_oligo, model_micro, model_neuron (each +1 single
  composition covariate). q-values no cutoff + leave-one-out. Single-covariate
  per model (cleaner than phase-1 multi-covariate; no over-parameterization).
- ART confound (Sarah's reframe): buckets C,D were on ART (HIV antiretroviral
  therapy), A,B,E were not. comp1 (2024neg-onART vs 2026neg-offART) characterizes
  ART/cohort/region, not toxo. comp2/comp3 are toxo contrasts confounded with
  ART (pos off-ART vs neg on-ART). comp4 (the "and the original" verbal addendum:
  pos vs 2026neg, both off-ART) is the patient-clean toxo contrast — trades the
  ART confound for a 2024-vs-2026 technical/cohort one. No bucket F (BG toxo-neg
  off-ART) can exist — grants are gone, no new data — so D's usability as a
  negative control for B is gauged by extrapolation from comp1, or BG is dropped.
- Result pattern: within-tissue toxo comparisons (2,3) — composition adjustment
  INCREASES significant tiles (acts as removable nuisance, sharpens toxo signal).
  comp4 (cleanest patient-axis) carries the most signal (1299 unadjusted),
  modest adjustment effect. Oligo cleanest covariate throughout; neuron most
  entangled (r=0.92 in comp1 → there it removes the contrast itself). LOO stable
  everywhere; 3-sample-arm comps (1,4) drop to min.per.group=2 on those drops.

---

## House style threaded through all of it

Data is holy; inputs reproducible/versioned; scripts self-contained and runnable
cold from disk (no session state — the scMD warm-session table is exactly the
anti-pattern). Loud failure over silent fallback (`ln -sn` not `-sf`; the
Guintivano gate stops rather than emitting a bad number; QC drives a settings
file that functions read else scream-and-die). Settings/manifests exported as
first-class deliverables shipped WITH results, in the same space as the tables.
Provenance written in phenomenon-language (TSS enrichment, percentMT, neuron
fraction), legible to a stranger in 2077 who's never heard of the tool. Route
around opacity rather than trust it ("pug metaphor": not angry at the dog, angry
at the breeding decision).
