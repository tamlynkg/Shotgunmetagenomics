Sed ‘s/>/>NK3_/g’ NK3scaffolds.fasta > ../binning/NK3scaffolds.fasta 

Cat *fasta > All_scaffolds.fa

Fasta_contig_size -s 1500 All_scaffolds.fa (Contig length:1500)

Mv All_scaffolds.fa All_scaffolds_1500+.fa 
