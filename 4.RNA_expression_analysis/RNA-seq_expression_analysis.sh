### remove the adapters and low quality reads

fastp -i ECI-2-120_C1_left.fastq -I ECI-2-120_C1_right.fastq -o ECI-2-120_C1_.clean.1.fastq -O ECI-2-120_C1_.clean.2.fastq -w 6


### Mapping the clean reads to the EC-S1 reference genome

hisat2 -x EC-S1_pseudomolecules_V1_1kb.fasta -1 ECI-2-120_C1_.clean.1.fastq -2 ECI-2-120_C1_.clean.2.fastq -p 6 |samtools view -@ 6 -Sb - | samtools sort -@ 6 -o ECI-2-120_C1_.sort.bam

### Count the mapped reads

Rscript run-featurecounts.R -b ECI-2-120_C1_.sort.bam -g new.gene.gtf -o ECI-2-120_C1

### Compare the RNA-seq expression pattern between samples

Rscript runDESeq2.R
