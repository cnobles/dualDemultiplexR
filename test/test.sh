Rscript dualDemultiplex.R -m test/sampleInfo.tsv \
  --read1 test/Data/Undetermined_S0_L001_R1_001.fastq.gz \
  --read2 test/Data/Undetermined_S0_L001_R2_001.fastq.gz \
  --index1 test/Data/Undetermined_S0_L001_I1_001.fastq.gz \
  --index2 test/Data/Undetermined_S0_L001_I2_001.fastq.gz \
  -o test/test_output --compress
