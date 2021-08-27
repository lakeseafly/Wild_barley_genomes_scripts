### Mapping the EC-N1 genome to EC-S1 genome and filter the mapping results 

~/biosoft/mummer-4.0.0beta2/nucmer -l 100 -c 500 -t 28 -p EC-S1_VS_EC-N1 EC-S1_pseudomolecules_V1_1kb_ref.fasta EC-N1_pseudomolecules_V1_1kb_ref.fasta
~/biosoft/mummer-4.0.0beta2/delta-filter -m -i 90 -l 100 EC-S1_VS_EC-N1.delta > EC-S1_VS_EC-N1.filter.delta
~/biosoft/mummer-4.0.0beta2/show-coords -THrd EC-S1_VS_EC-N1.filter.delta > EC-S1_VS_EC-N1.filter.coords

### Running Syri to detect the SV between EC-S1 and EC-N1 genome

python3.5 syri -c EC-S1_VS_EC-N1.filter.coords -d EC-S1_VS_EC-N1.filter.delta -r EC-S1_pseudomolecules_V1_1kb_ref.fasta -q EC-N1_pseudomolecules_V1_1kb_ref.fasta--nosnp --prefix EC-S1_VS_EC-N1.syri --lf EC-S1_VS_EC-N1_ref.log
