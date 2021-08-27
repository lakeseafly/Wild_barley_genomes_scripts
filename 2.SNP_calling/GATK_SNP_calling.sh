### Samples were mapped to the EC-S1 reference genome using BWA

bwa mem -M -R '@RG\tID:WB24\tLB:WB24\tPL:ILLUMINA\tPM:HISEQ\tSM:WB24' -t 24 EC-S1_pseudomolecules_V1_1kb.fasta WB24_1.fastq WB24_2.fastq | samtools view -@ 12 -Sb - |samtools sort -@ 12 -o WB24.sort.bam

### Finding and removing the duplicate reads in BAM file

sambamba markdup -r -t 16 WB24.sort.bam WB24.mark.bam

### HaplotypeCalling using GATK3.8

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R EC-S1_pseudomolecules_V1_1kb.fasta -I WB24.mark.bam -o WB24.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -nct 24 -allowPotentiallyMisencodedQuals

### After all the samples are processed the above steps, merge the generated g.vcf files of all the samples

java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R EC-S1_pseudomolecules_V1_1kb.fasta -V WD.gvcf.list -o WD_Bar.raw.vcf -newQual  -nt 28

### Extract the SNP and INDEL variants

java -jar ~/biosoft/GenomeAnalysisTK.jar -T SelectVariants -R EC-S1_pseudomolecules_V1_1kb.fasta -V WD_Bar.raw.vcf -selectType SNP -o SNPs.vcf

java -jar ~/biosoft/GenomeAnalysisTK.jar -T SelectVariants -R EC-S1_pseudomolecules_V1_1kb.fasta -V WD_Bar.raw.vcf -selectType INDEL  -o INDELs.vcf






