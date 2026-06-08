# For each of the 9 samples
# biscuit pileup -@ 16 /projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta CE167_GRCh38_no_alt_sorted.bam -o CE167_GRCh38.vcf
# biscuit vcf2bed -t cg CE167_GRCh38.vcf.gz > CE167_GRCh38_CpG.bed

for sample in CA192 CA346 CB239 CC249 CE167 CE234 LG30 LG31 LG52; do
    biscuit pileup -@ 16 -m 30 /projects1/references/GRCh38_no_alt_analysis_set/GRCh38_no_alt_analysis_set.fasta ${sample}/${sample}_GRCh38_no_alt_sorted_mq30.bam -o ${sample}/${sample}_GRCh38_mq30.vcf
    biscuit vcf2bed -t cg ${sample}/${sample}_GRCh38_mq30.vcf > ${sample}/${sample}_GRCh38_mq30_CpG.bed
done
