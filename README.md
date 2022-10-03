# Shotgunmetagenomics
Scripts used in the analysis of fasta files produced using Illumina

Analysis of shotgun metagenomics data

Part 1 ----

Quality control analysis of raw short reads using FASTQC tool kit and Prinseq tool.
Prinseq manual: http://prinseq.sourceforge.net/manual.html

Further based on the quality control report, data will be filtered by prinseq tool.
Removing the contamination from the sequencing data.

Download the PhiX and host genome from the NCBI database.

Decont.py script will be used to remove the expected contaminants from the sequencing data.

Part 2 ---

A)	BlastX analysis of short reads using “Diamond” tool and NCBI-NR database.
The Diamond output import into the MEGAN community edition (http://ab.inf.uni-tuebingen.de/data/software/megan6/download/manual.pdf) for taxonomic and functional analyses.

B)	Short reads assembly using metaSPAdes assembler
Quality check of the assembly by MetaQuast pipeline.
Best contigs will be further used for the ORF prediction using Prodigal pipelines.
Blastp analysis of ORFs using NCBI blastp tool and NCBI-NR protein database.
Analysis of functional ORFs in MEGAN CE edition.

Part 3 ----

Statistical comparison of the dataset obtained from the different sample types
This analysis can be done in the STAMP tool
Download STAMP and read the manual: http://kiwi.cs.dal.ca/Software/STAMP

Genome reconstruction from metagenomes:
Assembled contigs will be used for the genome reconstruction.
Reads will be aligned to the contigs (longer than 1500 bp) using bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) aligner. The alignment will be sorted as a bam file using samtools (http://www.htslib.org/doc/samtools-1.2.html).
Binning tools will be used for the de novo binning of the genomes.
Metabat (https://bitbucket.org/berkeleylab/metabat),
Maxbin (https://downloads.jbei.org/data/microbial_communities/MaxBin/README.txt),
MyCC (http://download2.nust.na/pub4/sourceforge/s/sb/sb2nhri/MyCC/manual%20of%20MyCC.pdf)
Refinement of the genome bins
The reconstructed genomes will used for the quality assessment using checkm tool (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/pdf/1043.pdf)
Further the multiple binnings will be combined for the refinement process (remove contamination, duplicate sequences and completeness) using binning refiner (https://github.com/songweizhi/Binning_refiner/blob/master/README.md) and DAS tools (https://github.com/cmks/DAS_Tool).
Quality of the reconstructed genomes will be assessd at each step using checkm tool.
