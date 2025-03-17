Scripts needed to extract relavant data from BAM files and format for MLM. 

1. Extract all reads that span position of variant in each of Bam files of interest with `samtools view`
2. Run Generate_Fragment_Files.py on the Sam/Bam files to identify the mutant and wild type fragments and collect relevant features.
3. Run CHIP_fragmentLevel_sumStats.R to calculate summary statistics on the fragments. 
