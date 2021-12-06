### Running EDTA to identify the EC-S1 and EC-N1  de novo transposable element (TE) libraries 

EDTA.pl --genome fixorder_Damai_S_pseudomolecules_V1.fasta --anno 1 --evaluate 1 --threads 28
EDTA.pl --genome fixorder_Damai_N_pseudomolecules_V1.fasta --anno 1 --evaluate 1 --threads 28

### Remove the highly similiar BARE sequences
cd-hit -i barley.BARE.TE.fasta -o barley.BARE.TE.non-redundant.fasta -c 0.9

### Using BLAST+ to search for the EC-S1 and EC-N1 specific BARE sequences and filtered the results to generate the final 

blastn --db barley.BARE.TE.non-redundant.fasta --query EC-S1.TE.fa --out EC-S1.BARE.results.out --outfmt 6 --evalue 0.00001
blastn --db barley.BARE.TE.non-redundant.fasta --query EC-N1.TE.fa --out EC-N1.BARE.results.out --outfmt 6 --evalue 0.00001


