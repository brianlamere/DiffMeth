#!/bin/bash
# ToxoDM alignment — option 1: GRCh38_no_alt_analysis_set, no HXB2 (matches project 1)
# Aligns the 14 NEW samples (10 BG + 4 FMC); symlinks the 9 inherited cortex/LG beds.
# Settings copied exactly from project-1 @PG: biscuit 1.8.0, default MQ40 pileup, vcf2bed -t cg.
set -euo pipefail

REF=/projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta
FQ=/projects1/ToxoDM/fastqs
OUT=/projects1/ToxoDM/working/align
ATHREADS=16   # align threads (matches project 1 -@ 16)
STHREADS=14   # sort threads (matches project 1 -@ 14)

align_one () {   # align_one <cohort> <tissue> <sample>
    local cohort=$1 tissue=$2 sample=$3
    local sdir="$FQ/$cohort/$tissue/$sample"
    local odir="$OUT/$cohort/$tissue/$sample"
    mkdir -p "$odir"
    # R1/R2 are symlinks named R1.* / R2.* (.gz or .fq); glob picks them up
    local r1=$(echo "$sdir"/R1.*) r2=$(echo "$sdir"/R2.*)
    local tag="${sample}_${tissue}"
    echo ">> aligning $cohort/$tissue/$sample"
    biscuit align -@ $ATHREADS "$REF" "$r1" "$r2" > "$odir/${tag}_raw.bam"
    samtools sort -@ $STHREADS "$odir/${tag}_raw.bam" > "$odir/${tag}_sorted.bam"
    samtools index "$odir/${tag}_sorted.bam"
    biscuit pileup -@ $ATHREADS "$REF" "$odir/${tag}_sorted.bam" -o "$odir/${tag}.vcf"
    biscuit vcf2bed -t cg "$odir/${tag}.vcf" > "$odir/${tag}_GRCh38_CpG.bed"
    rm -f "$odir/${tag}_raw.bam"   # raw unsorted bam is disposable once sorted exists
    echo "   done: $odir/${tag}_GRCh38_CpG.bed"
}

echo "=============================================="
echo " NEW ALIGNMENTS (14): 2024 BG (10) + FMC (4)"
echo "=============================================="
# 2024 BG — 6 positives + 4 negatives
for s in CA192 CA346 CB239 CC249 CE167 CE234 Atpsy_15 Atpsy_20 Atpsy_21 LG29; do
    align_one 2024 BG "$s"
done
# 2024 FMC — 4 negatives
for s in Atpsy_15 Atpsy_20 Atpsy_21 LG29; do
    align_one 2024 FMC "$s"
done

echo "=============================================="
echo " INHERITED BEDS (9): symlink project-1 CC beds"
echo "=============================================="
# 6 CC positives from 240131, 3 LG from 260202 — already GRCh38_no_alt, biscuit 1.8.0
declare -A CC_SRC=(
  [CA192]=/projects2/240131_LH00444_0048_A22GLTCLT3/CA192_CC/CA192_GRCh38_CpG.bed
  [CA346]=/projects2/240131_LH00444_0048_A22GLTCLT3/CA346_CC/CA346_GRCh38_CpG.bed
  [CB239]=/projects2/240131_LH00444_0048_A22GLTCLT3/CB239_CC/CB239_GRCh38_CpG.bed
  [CC249]=/projects2/240131_LH00444_0048_A22GLTCLT3/CC249_CC/CC249_GRCh38_CpG.bed
  [CE167]=/projects2/240131_LH00444_0048_A22GLTCLT3/CE167_CC/CE167_GRCh38_CpG.bed
  [CE234]=/projects2/240131_LH00444_0048_A22GLTCLT3/CE234_CC/CE234_GRCh38_CpG.bed
)
for s in "${!CC_SRC[@]}"; do
    odir="$OUT/2024/CC/$s"; mkdir -p "$odir"
    ln -sn "${CC_SRC[$s]}" "$odir/${s}_CC_GRCh38_CpG.bed"
    echo "   linked 2024/CC/$s -> ${CC_SRC[$s]}"
done
declare -A LG_SRC=(
  [LG30]=/projects2/260202_LH00444_0467_A233LHLLT3/LG30/LG30_GRCh38_CpG.bed
  [LG31]=/projects2/260202_LH00444_0467_A233LHLLT3/LG31/LG31_GRCh38_CpG.bed
  [LG52]=/projects2/260202_LH00444_0467_A233LHLLT3/LG52/LG52_GRCh38_CpG.bed
)
for s in "${!LG_SRC[@]}"; do
    odir="$OUT/2026/CC/$s"; mkdir -p "$odir"
    ln -sn "${LG_SRC[$s]}" "$odir/${s}_CC_GRCh38_CpG.bed"
    echo "   linked 2026/CC/$s -> ${LG_SRC[$s]}"
done

echo ""
echo "All done. New beds + inherited beds under $OUT"
