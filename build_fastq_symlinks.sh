#!/bin/bash
# Build /projects1/ToxoDM/fastqs/<cohort>/<tissue>/<sample>/{R1,R2} symlink tree
# Sources: toxo3 val_bsfiles (2024 samples) + 260202 trimmed (2026 LG samples)
# Clean R1/R2 names regardless of messy source filenames.
set -euo pipefail

BASE=/projects1/ToxoDM/fastqs
TOXO3=/projects2/toxo3/val_bsfiles
LG2026=/projects2/260202_LH00444_0467_A233LHLLT3

mklink () {  # mklink <cohort> <tissue> <sample> <r1_src> <r2_src>
    local d="$BASE/$1/$2/$3"
    mkdir -p "$d"
    #ln -sfn "$4" "$d/R1.${4##*.}"
    ln -sn "$4" "$d/R1.${4##*.}"
    #ln -sfn "$5" "$d/R2.${5##*.}"
    ln -sn "$5" "$d/R2.${5##*.}"
    echo "  $1/$2/$3 -> $(basename "$4"), $(basename "$5")"
}

echo "== 2024 cohort: toxo3 val_bsfiles =="
# 2024 positives (toxo+): CC and BG for the 6 distinct donors
# NOTE: source dirs use CC_249 / CE_167 (underscored) in CC tissue, CC249/CE167 in BG
# 2024 CC positives
mklink 2024 CC CA192  "$TOXO3/CC/CA192/CA192_Cerebral_Cortex_S14_L007_R1_001_val_1.fq.gz"  "$TOXO3/CC/CA192/CA192_Cerebral_Cortex_S14_L007_R2_001_val_2.fq.gz"
mklink 2024 CC CA346  "$TOXO3/CC/CA346/CA346_Cerebral_Cortex_S15_L007_R1_001_val_1.fq.gz"  "$TOXO3/CC/CA346/CA346_Cerebral_Cortex_S15_L007_R2_001_val_2.fq.gz"
mklink 2024 CC CB239  "$TOXO3/CC/CB239/CB239_Cerebral_Cortex_S18_L007_R1_001_val_1.fq.gz"  "$TOXO3/CC/CB239/CB239_Cerebral_Cortex_S18_L007_R2_001_val_2.fq.gz"
mklink 2024 CC CC249  "$TOXO3/CC/CC_249/CC_249_Cerebral_Cortex_S16_L007_R1_001_val_1.fq.gz" "$TOXO3/CC/CC_249/CC_249_Cerebral_Cortex_S16_L007_R2_001_val_2.fq.gz"
mklink 2024 CC CE167  "$TOXO3/CC/CE_167/CE_167_Cerebral_Cortex_S17_L007_R1_001_val_1.fq.gz" "$TOXO3/CC/CE_167/CE_167_Cerebral_Cortex_S17_L007_R2_001_val_2.fq.gz"
mklink 2024 CC CE234  "$TOXO3/CC/CE234/CE234_Cerebral_Cortex_S8_L007_R1_001_val_1.fq.gz"   "$TOXO3/CC/CE234/CE234_Cerebral_Cortex_S8_L007_R2_001_val_2.fq.gz"
# 2024 BG positives
mklink 2024 BG CA192  "$TOXO3/BG/CA192/CA192_Basal_Ganglia_S7_L007_R1_001_val_1.fq.gz"  "$TOXO3/BG/CA192/CA192_Basal_Ganglia_S7_L007_R2_001_val_2.fq.gz"
mklink 2024 BG CA346  "$TOXO3/BG/CA346/CA346_Basal_Ganglia_S1_L007_R1_001_val_1.fq.gz"  "$TOXO3/BG/CA346/CA346_Basal_Ganglia_S1_L007_R2_001_val_2.fq.gz"
mklink 2024 BG CB239  "$TOXO3/BG/CB239/CB239_Basal_Ganglia_S4_L007_R1_001_val_1.fq.gz"  "$TOXO3/BG/CB239/CB239_Basal_Ganglia_S4_L007_R2_001_val_2.fq.gz"
mklink 2024 BG CC249  "$TOXO3/BG/CC249/CC249_Basal_Ganglia_S2_L007_R1_001_val_1.fq.gz"  "$TOXO3/BG/CC249/CC249_Basal_Ganglia_S2_L007_R2_001_val_2.fq.gz"
mklink 2024 BG CE167  "$TOXO3/BG/CE167/CE167_Basal_Ganglia_S3_L007_R1_001_val_1.fq.gz"  "$TOXO3/BG/CE167/CE167_Basal_Ganglia_S3_L007_R2_001_val_2.fq.gz"
mklink 2024 BG CE234  "$TOXO3/BG/CE234/CE234_Basal_Ganglia_S20_L007_R1_001_val_1.fq.gz" "$TOXO3/BG/CE234/CE234_Basal_Ganglia_S20_L007_R2_001_val_2.fq.gz"
# 2024 negatives (toxo-): FMC
mklink 2024 FMC Atpsy_15 "$TOXO3/FMC/Atpsy_15/Atpsy_15_Frontal_Motor_Cortex_S23_L007_R1_001_val_1.fq.gz" "$TOXO3/FMC/Atpsy_15/Atpsy_15_Frontal_Motor_Cortex_S23_L007_R2_001_val_2.fq.gz"
mklink 2024 FMC Atpsy_20 "$TOXO3/FMC/Atpsy_20/Atpsy_20_Frontal_Motor_Cortex_S22_L007_R1_001_val_1.fq.gz" "$TOXO3/FMC/Atpsy_20/Atpsy_20_Frontal_Motor_Cortex_S22_L007_R2_001_val_2.fq.gz"
mklink 2024 FMC Atpsy_21 "$TOXO3/FMC/Atpsy_21/Atpsy_21_Frontal_Motor_Cortex_S21_L007_R1_001_val_1.fq.gz" "$TOXO3/FMC/Atpsy_21/Atpsy_21_Frontal_Motor_Cortex_S21_L007_R2_001_val_2.fq.gz"
mklink 2024 FMC LG29     "$TOXO3/FMC/LG29/LG29_Frontal_Motor_Cortex_S24_L007_R1_001_val_1.fq.gz"         "$TOXO3/FMC/LG29/LG29_Frontal_Motor_Cortex_S24_L007_R2_001_val_2.fq.gz"
# Also BG versions of the 2024 negatives (comparison 3 needs neg BG)
mklink 2024 BG Atpsy_15 "$TOXO3/BG/Atpsy_15/Atpsy_15_Basal_Ganglia_S27_L007_R1_001_val_1.fq.gz" "$TOXO3/BG/Atpsy_15/Atpsy_15_Basal_Ganglia_S27_L007_R2_001_val_2.fq.gz"
mklink 2024 BG Atpsy_20 "$TOXO3/BG/Atpsy_20/Atpsy_20_Basal_Ganglia_S25_L007_R1_001_val_1.fq.gz" "$TOXO3/BG/Atpsy_20/Atpsy_20_Basal_Ganglia_S25_L007_R2_001_val_2.fq.gz"
mklink 2024 BG Atpsy_21 "$TOXO3/BG/Atpsy_21/Atpsy_21_Basal_Ganglia_S26_L007_R1_001_val_1.fq.gz" "$TOXO3/BG/Atpsy_21/Atpsy_21_Basal_Ganglia_S26_L007_R2_001_val_2.fq.gz"
mklink 2024 BG LG29     "$TOXO3/BG/LG29/LG29_Basal_Ganglia_S28_L007_R1_001_val_1.fq.gz"         "$TOXO3/BG/LG29/LG29_Basal_Ganglia_S28_L007_R2_001_val_2.fq.gz"

echo "== 2026 cohort: 260202 trimmed (negatives, CC) =="
mklink 2026 CC LG30 "$LG2026/LG30/trimmed/LG30_RRBS_S1_L003_R1_001_val_1.fq" "$LG2026/LG30/trimmed/LG30_RRBS_S1_L003_R2_001_val_2.fq"
mklink 2026 CC LG31 "$LG2026/LG31/trimmed/LG31_RRBS_S2_L003_R1_001_val_1.fq" "$LG2026/LG31/trimmed/LG31_RRBS_S2_L003_R2_001_val_2.fq"
mklink 2026 CC LG52 "$LG2026/LG52/trimmed/LG52_RRBS_S3_L003_R1_001_val_1.fq" "$LG2026/LG52/trimmed/LG52_RRBS_S3_L003_R2_001_val_2.fq"

echo ""
echo "Done. Tree:"
find "$BASE" -type l | sort
