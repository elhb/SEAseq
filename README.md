##Readme for the SEAseq scripts collection
This is a set of script for analysis of SEAseq data in fastq files.

###Example usage:

Initiate the analysis folder:
`script/SEAseq init -path tasks/testing`

Add some fastq files to be analyzed:
`script/SEAseq addfqs -path tasks/testing/ -r1 data/20130205_BC1DF7ACXX_HiSeq2000Spikein/P412_102_index12/130205_BC1DF7ACXX/4_130205_BC1DF7ACXX_P412_102_index12_1.fastq -r2 data/20130205_BC1DF7ACXX_HiSeq2000Spikein/P412_102_index12/130205_BC1DF7ACXX/4_130205_BC1DF7ACXX_P412_102_index12_2.fastq` 
`script/SEAseq addfqs -path tasks/testing/ -r1 data/20130205_BC1DF7ACXX_HiSeq2000Spikein/P412_101_index10/130205_BC1DF7ACXX/4_130205_BC1DF7ACXX_P412_101_index10_1.fastq -r2 data/20130205_BC1DF7ACXX_HiSeq2000Spikein/P412_101_index10/130205_BC1DF7ACXX/4_130205_BC1DF7ACXX_P412_101_index10_2.fastq`

Cluster the barcode sequences in the reads files:
`script/SEAseq clusterbarcodes -path tasks/testing/ -hm 4 -bm 4 -p 8`

Sort the reads to clusters:
`script/SEAseq sortreads -path tasks/testing/ -p 8 -stop 14869300`
