Mk metabat

#create bam file
Samtools view – bS MK2.sam > MK2.bam 

#The alignment was sorted as a bam file using samtools
#sort bam files
Samtools sort MK2.bam -o MK2.bam-sorted.bam 
