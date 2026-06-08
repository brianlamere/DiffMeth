# For each of the 9 samples
for sample in CA192 CA346 CB239 CC249 CE167 CE234 LG30 LG31 LG52; do
    samtools view -q 30 -@ 8 -b \
        ${sample}/${sample}_GRCh38_no_alt_sorted.bam \
        | samtools sort -@ 8 \
        -o ${sample}/${sample}_GRCh38_no_alt_sorted_mq30.bam
    samtools index ${sample}/${sample}_GRCh38_no_alt_sorted_mq30.bam
done
