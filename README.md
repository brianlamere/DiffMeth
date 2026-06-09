# DiffMeth

bucket 1 was built as:
 trim_galore --rrbs --paired --cores 14 -o CA192_CC/trimmed CA192_CC/CA192_Cerebral_Cortex_S14_L007_R1_001.fastq.gz CA192_CC/CA192_Cerebral_Cortex_S14_L007_R2_001.fastq.gz
 biscuit pileup -@ 16 /projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta CA192_CC/CA192_GRCh38_no_alt_sorted.bam -o CA192_CC/CA192_GRCh38.vcf.gz
 samtools sort -@ 12 CA192_CC/CA192_GRCh38-no-alt_raw.bam > CA192_CC/CA192_GRCh38-no-alt_sorted.bam


bucket 2 (2026) was built as:

 biscuit align -@ 16 /projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta LG30/trimmed/LG30_RRBS_S1_L003_R1_001_val_1.fq LG30/trimmed/LG30_RRBS_S1_L003_R2_001_val_2.fq > LG30/LG30_GRCh38_no_alt_raw.bam
 samtools sort -@ 14 LG30/LG30_GRCh38_no_alt_raw.bam > LG30/LG30_GRCh38_no_alt_sorted.bam
