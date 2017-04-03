#!/bin/R
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("yaml"))
suppressMessages(library("pander"))

code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

#' Set up and gather command line arguments
parser <- ArgumentParser(description = "R-based demultiplexing for GuideSeq.")
parser$add_argument(
  "-m", "--manifest", type = "character", help = "GuideSeq manifest yaml.")
parser$add_argument(
  "--read1", type = "character", help = "Path to Illumina R1 file (FASTQ).")
parser$add_argument(
  "--read2", type = "character", help = "Path to Illumina R2 file (FASTQ).")
parser$add_argument(
  "--index1", type = "character", help = "Path to Illumina I1 file (FASTQ).")
parser$add_argument(
  "--index2", type = "character", help = "Path to Illumina I2 file (FASTQ).")
parser$add_argument(
  "-o", "--outfolder", nargs = 1, type = "character", help = "Output folder.")
parser$add_argument(
  "-p", "--poolreps", action = "store_true", help = "Pools replicates.")
parser$add_argument(
  "--maxMismatch", nargs = 1, type = "integer", help = "Max mismatch allowed in barcodes.")
parser$add_argument(
  "--barcode1Length", nargs = 1, type = "integer", default = 8, 
  help = "Length of barcode1, in nucleotides. Default = 8.")
parser$add_argument(
  "--barcode2Length", nargs = 1, type = "integer", default = 8,
  help = "Length of barcode2, in nucleotides. Default = 8.")
parser$add_argument(
  "-c", "--cores", nargs = 1, type = "integer", help = "Max cores to be used.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("manifest :", "index1 :", "index2 :", "read1 :", "read2 :", 
          "outfolder :", "poolreps :", "cores :", "maxMismatch :", 
          "barcode1Length :", "barcode2Length :"),
        input_table$Variables),]
pandoc.title("Demultiplex Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Load additional dependencies
add_dependencies <- c("stringr", "ShortRead", "Biostrings")

null <- suppressMessages(
  sapply(add_dependencies, library, character.only = TRUE))

rm(null)

# Load guideseq manifest file
manifest <- yaml.load_file(args$manifest)
samples.df <- data.frame(
  "sample" = names(manifest$samples),
  "description" = sapply(manifest$samples, function(x) x$description),
  "barcode1" = sapply(manifest$samples, function(x) x$barcode1),
  "barcode2" = sapply(manifest$samples, function(x) x$barcode2),
  row.names = NULL
)

# Read in I1 and I2 sequences
parseIndexReads <- function(barcode, indexFilePath, barcodeLength, maxMismatch){
  dependencies <- c("stringr", "ShortRead", "Biostrings")
  null <- sapply(dependencies, library, character.only = TRUE)
  index <- readFastq(indexFilePath)
  indexReads <- DNAStringSet(sread(index), start = 1, end = barcodeLength)
  names(indexReads) <- sapply(strsplit(as.character(id(index)), " "), "[[", 1)
  mindex <- vmatchPattern(
    barcode, 
    indexReads, 
    max.mismatch = maxMismatch)
  names(unlist(mindex))
}

cluster <- makeCluster(min(c(detectCores(), args$cores)))

I1.parsed <-  parLapply(
  cluster,
  unique(samples.df$barcode1), 
  parseIndexReads,
  indexFilePath = args$index1,
  barcodeLength = args$barcode1Length,
  maxMismatch = args$maxMismatch)

I2.parsed <- parLapply(
  cluster, 
  unique(samples.df$barcode2), 
  parseIndexReads,
  indexFilePath = args$index2,
  barcodeLength = args$barcode2Length,
  maxMismatch = args$maxMismatch)

stopCluster(cluster)

names(I1.parsed) <- unique(samples.df$barcode1)
names(I2.parsed) <- unique(samples.df$barcode2)

pandoc.title("Barcode1 breakdown:")
pandoc.table(sapply(I1.parsed, length), split.tables = Inf)
pandoc.title("Barcode2 breakdown:")
pandoc.table(sapply(I2.parsed, length), split.tables = Inf)

demultiplexedReadNames <- mapply(
  function(barcode1, barcode2){
    intersect(I1.parsed[[barcode1]], I2.parsed[[barcode2]])
  },
  barcode1 = samples.df$barcode1,
  barcode2 = samples.df$barcode2)
names(demultiplexedReadNames) <- paste0(
  samples.df$barcode1, samples.df$barcode2)

# As there is some flexibility in the barcode matching, some reads may be 
# be assigned to multiple samples. These reads are ambiguous and will be 
# removed.
ambiguousReads <- unique(
  unlist(demultiplexedReadNames)[duplicated(unlist(demultiplexedReadNames))])
demultiplexedReadNames <- lapply(demultiplexedReadNames, function(x, reads){
  x[!x %in% reads]},
  reads = ambiguousReads)

# Reads by sample
samples.df$read_counts <- sapply(demultiplexedReadNames, length)
pandoc.title("Read counts for each sample")
pandoc.table(samples.df, split.tables = Inf)
# Ambiguous reads
message(paste0("Ambiguous reads: ", length(ambiguousReads)))
# Unassigned reads
readNames <- id(readFastq(args$index1))
readNames <- sapply(strsplit(as.character(readNames), " "), "[[", 1)
unassignedReadNames <- readNames[
  !readNames %in% unlist(demultiplexedReadNames)]
unassignedReadNames <- unassignedReadNames[
  !unassignedReadNames %in% ambiguousReads]
message(paste0("\nUnassigned reads: ", length(unassignedReadNames)))

# Create multiplex dataframe for subseting sequencing files
multiplexedData <- data.frame(
  "sample" = Rle(
    values = samples.df$sample, 
    length = sapply(demultiplexedReadNames, length)),
  "readName" = unlist(demultiplexedReadNames),
  row.names = NULL)

ambiguousData <- data.frame(
  "sample" = rep("ambiguous", length(ambiguousReads)),
  "readName" = ambiguousReads,
  row.names = NULL)

unassignedData <- data.frame(
  "sample" = rep("unassigned", length(unassignedReadNames)),
  "readName" = unassignedReadNames,
  row.names = NULL)

multiplexedData <- rbind(multiplexedData, ambiguousData, unassignedData)

stopifnot(nrow(multiplexedData) == length(readNames))

if(args$poolreps){
  multiplexedData$sample <- sapply(
    strsplit(multiplexedData$sample, "-"), "[[", 1)
}

message(paste0("\nReads writen to files: ", nrow(multiplexedData)))

# Write files to read files to outfolder directory
writeDemultiplexedSequences <- function(multiplexedData, readFilePath, 
                                        type, outfolder){
  reads <- readFastq(readFilePath)
  ids <- sapply(strsplit(as.character(id(reads)), " "), "[[", 1)
  reads <- reads[match(multiplexedData$readName, ids)]
  reads <- split(reads, multiplexedData$sample)
  
  null <- lapply(1:length(reads), function(i, reads, type, outfolder){
    filePath <- file.path(
      outfolder, paste0(names(reads[i]), ".", type, ".fastq"))
    writeFastq(reads[[i]], file = filePath)
    message(paste0("\nWrote ", length(reads[[i]]), " reads to:\n", filePath, "."))
  }, reads = reads, type = type, outfolder = outfolder)
  return(list(readFilePath, type, outfolder))
}

demultiplex.I1 <- writeDemultiplexedSequences(
  multiplexedData, args$index1, "i1", args$outfolder)
demultiplex.I2 <- writeDemultiplexedSequences(
  multiplexedData, args$index2, "i2", args$outfolder)
demultiplex.R1 <- writeDemultiplexedSequences(
  multiplexedData, args$read1, "r1", args$outfolder)
demultiplex.R2 <- writeDemultiplexedSequences(
  multiplexedData, args$read2, "r2", args$outfolder)

message("\nDemultiplexing complete.")
q()
