# dualDemultiplexR
Demultiplex Illumina FASTQ files using dual barcodes (I1 and I2).

Usage:
```
Rscript path/to/dualDemultiplex.R -m manifest.yml \
  --read1 READ1.fastq --read2 READ2.fastq --index1 INDEX1.fastq --index2 INDEX2.fastq \
  --outFolder ~/demultiplexed

Rscript path/to/dualDemultiplex.R -m manifest.yml \
  --read1 READ1.fastq --read2 READ2.fastq --index1 INDEX1.fastq --index2 INDEX2.fastq \
  --outFolder ~/demultiplexed --poolreps --maxMismatch 1 --barcode1Length 8 --barcode2Length 8 \
  --readNamePattern [\\w:-]+ --compress --cores 4
```

## Sample manifest format
Manifests contain the barcode information for each sample, currently in a yaml format. An example is show below:
```
# This is an example sample manifest.
samples :
    treated_sample-1 :
        barcode1 : TAAGGCGA
        barcode2 : TAGATCGC
        description : Replicate 1 of treated_sample
    
    treated_sample-2 :
        barcode1 : TAAGGCGA
        barcode2 : CGATCCTA
        description : Replicate 2 of treated_sample
        
    alt_treated_sample-1 :
        barcode1 : NATCGTCA
        barcode2 : NCTGGTAC
        description : Replicate 1 of alt_treated_sample
        
...
```
If replicates are included and are to be pooled, sample names should be in the above format (sampleName-#). Pooled replicates with be consolidated if the "-p" or "--poolreps" flag is used with the "-" as the delimiter between sampleNames and replicate designations. Likewise, sampleNames should not include the "-" symbol if pooling replicates is desired.

Barcode sequences may contain ambiguous nucleotides. It is recommended that the user chooses to use either ambiguous nucleotides or declare the maximal number of mismatches allowed. Using both at the same time may lead to undesired demultiplexing. Nucleotide substitution matrix NUC4.4 used for scoring (ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4). Both ambiguous nucleotide matching and maxMismatch allowance can be useful in their own way. For poor sequencing quality, ambigous nucleotide demultiplexing may be able to acquire more reads that contained low quality index reads, while maxMismatch may help if phasing/pre-phasing is a major issue in the sequencing file.

## Arguments
**[-h, --help]** Help information regarding input format and arguments available.

**[--read1, --read2, --index1, --index2]** File paths to the respective FASTQ files.

**[-o, --outFolder]** Output directory for demultiplexed files.

**[-p, --poolreps]** Pools replicates. Replicate designation delimited by "-". ie. sample-1, sample-2, ...

**[--maxMismatch]** Allowable mismatch in barcode sequences. 

**[--barcode1Length, --barcode2Length]** Length of barcode sequences for demultiplexing. Default is 8 nt.

**[--readNamePattern]** Regex pattern to capture read names without read-type specific info.

**[--compress]** Output fastq files are gzipped.

**[-c, --cores]** Number of maximum cores to parallel the processing during certain steps.

## Dependencies
dualDemultiplexR is coded in R, and was developed on v3.2.2, though it should run with earlier versions given the appropriate dependencies. The script uses 6 additional packages:
  * argparse
  * yaml
  * pander
  * ShortRead
  * Biostrings
  * stringr
  * parallel (if multicore processing is desired)
