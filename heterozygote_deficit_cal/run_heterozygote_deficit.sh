### using the plink to calculate the Hobs, Hexp and F value
plink --vcf SNPs.vcf --hardy 

### using the windowscanr to estimate the F value by 100kb windows

Rscript Hobs_sliding.r SNP_plink.hwe 100000 10000
