# DiffMeth

Start with deconv_buckets.R, then use annotate.R, then use diffmeth.R

buckets filled as:
``` 
 trim_galore --rrbs --paired --cores 14 -o CA192_CC/trimmed CA192_CC/CA192_Cerebral_Cortex_S14_L007_R1_001.fastq.gz CA192_CC/CA192_Cerebral_Cortex_S14_L007_R2_001.fastq.gz
 biscuit align -@ 16 /projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta CA192_CC/trimmed/CA192_Cerebral_Cortex_S14_L007_R1_001_val_1.fq CA192_CC/trimmed/CA192_Cerebral_Cortex_S14_L007_R2_001_val_2.fq > CA192_CC/CA192_GRCh38-no-alt_raw.bam
 biscuit pileup -@ 16 /projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta CA192_CC/CA192_GRCh38_no_alt_sorted.bam -o CA192_CC/CA192_GRCh38.vcf.gz
 biscuit vcf2bed -t cg CA192_CC/CA192_GRCh38.vcf.gz > CA192_CC/CA192_GRCh38_CpG.bed
 samtools sort -@ 12 CA192_CC/CA192_GRCh38-no-alt_raw.bam > CA192_CC/CA192_GRCh38-no-alt_sorted.bam
```

