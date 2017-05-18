#!/bin/R
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))

code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

#' Set up and gather command line arguments
parser <- ArgumentParser(
  description = "R-based demultiplexing for dual-barcoded Illumina sequencing runs.")
parser$add_argument(
  "-m", "--manifest", type = "character", help = "Manifest file (*.yaml, *.yml, *.csv, *.tsv).")
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
  "--singleBarcode", nargs = 1, type = "character", 
  help = "Demultiplex with only a single barcode. Manifest should contain empty data for unused barcode. (Default = FALSE, options include 'index1' and 'index2'.)")
parser$add_argument(
  "--maxMismatch", nargs = 1, type = "integer", default = 0,
  help = "Max mismatch allowed in barcodes (Default = 0). It is recommended to use either ambiguous nucleotide codes or maxMismatch, not both.")
parser$add_argument(
  "--barcode1Length", nargs = 1, type = "integer", default = 8, 
  help = "Length of barcode1, in nucleotides. Default = 8.")
parser$add_argument(
  "--barcode2Length", nargs = 1, type = "integer", default = 8,
  help = "Length of barcode2, in nucleotides. Default = 8.")
parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", default = "[\\w:-]+",
  help = "Regular expression for pattern matching read names. Should not contain R1/R2/I1/I2 specific components. Default is [\\w:-]+")
parser$add_argument(
  "--compress", action = "store_true", help = "Output fastq files are gzipped.")
parser$add_argument(
  "-c", "--cores", nargs = 1, default = 0, type = "integer", 
  help = "Max cores to be used. If 0 (default), program will not utilize parallel processing.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if(!args$singleBarcode %in% c(FALSE, "index1", "index2")){
  stop("Please choose either 'index1' or 'index2' for single demultiplexing.")
}

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("manifest :", "index1 :", "index2 :", "read1 :", "read2 :", 
          "outfolder :", "poolreps :", "singleBarcode :", "cores :", "maxMismatch :", 
          "barcode1Length :", "barcode2Length :", "readNamePattern :"),
        input_table$Variables),]
pandoc.title("Demultiplex Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Create output directory if not currently available
if(!file.exists(args$outfolder)){
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.")
}

# Load additional dependencies
add_dependencies <- c("stringr", "ShortRead", "Biostrings")
addDependsLoaded <- suppressMessages(
  sapply(add_dependencies, require, character.only = TRUE))
if(!all(addDependsLoaded)) stop("Check dependancies: 'stringr', 'Biostrings', 'ShortRead'.")
source(file.path(
  code_dir, "supporting_scripts/binary_ambiguous_nucleotide_scoring_matrix.R"))

submat <- banmat()

# Load guideseq manifest file
fileExt <- unlist(strsplit(args$manifest, "\\."))
fileExt <- fileExt[length(fileExt)]

if(fileExt %in% c("yaml", "yml")){
  suppressMessages(library("yaml"))
  if(!"package:yaml" %in% search()) stop("Package:yaml not loaded or installed.")
  manifest <- yaml.load_file(args$manifest)
  samples.df <- data.frame(
    "sampleName" = names(manifest$samples),
    "barcode1" = sapply(manifest$samples, function(x) x$barcode1),
    "barcode2" = sapply(manifest$samples, function(x) x$barcode2),
    row.names = NULL)
}else if(fileExt == "csv"){
  manifest <- read.csv(args$manifest)
  samples.df <- manifest[, c("sampleName", "barcode1", "barcode2")]
}else if(fileExt == "tsv"){
  manifest <- read.delim(args$manifest)
  samples.df <- manifest[, c("sampleName", "barcode1", "barcode2")]
}

if(args$singleBarcode == FALSE){
  uniqueSamples <- nrow(samples.df[,c("barcode1", "barcode2")]) == 
    nrow(unique(samples.df[,c("barcode1", "barcode2")]))
  if(!uniqueSamples) stop("Samples do not have unique barcodes. Please correct.")
}else{
  if(args$singleBarcode == "index1"){
    uniqueSamples <- length(samples.df[,c("barcode1")]) == 
      length(unique(samples.df[,"barcode1"]))
    if(!uniqueSamples) stop("Samples do not have unique barcodes. Please correct.")
  }else{
    uniqueSamples <- length(samples.df[,c("barcode2")]) == 
      length(unique(samples.df[,"barcode2"]))
    if(!uniqueSamples) stop("Samples do not have unique barcodes. Please correct.")
  }
}

# Read in I1 and I2 sequences
parseIndexReads <- function(barcode, indexFilePath, barcodeLength, maxMismatch, 
                            submat, readNamePattern){
  # Require packages for parallel processing
  dependencies <- c("stringr", "ShortRead", "Biostrings")
  loaded <- sapply(dependencies, require, character.only = TRUE)
  stopifnot(all(loaded))

  # Load index file sequences and sequence names
  index <- readFastq(indexFilePath)
  indexReads <- DNAStringSet(sread(index), start = 1, end = barcodeLength)
  names(indexReads) <- sapply(strsplit(as.character(id(index)), " "), "[[", 1)

  # Identify read names with sequences above or equal to the minscore
  scores <- pairwiseAlignment(
    indexReads, barcode, type = "global", substitutionMatrix = submat, 
    gapOpening = 3, gapExtension = 1, scoreOnly = TRUE)
  str_extract(as.character(id(index)), readNamePattern)[
    scores >= (barcodeLength - maxMismatch)]
}

if(args$cores > 0){
  suppressMessages(library(parallel))
  cluster <- makeCluster(min(c(detectCores(), args$cores)))
  
  if(args$singleBarcode == FALSE | args$singleBarcode == "index1"){
    I1.parsed <-  parLapply(
      cluster,
      unique(samples.df$barcode1), 
      parseIndexReads,
      indexFilePath = args$index1,
      barcodeLength = args$barcode1Length,
      maxMismatch = args$maxMismatch,
      submat = submat,
      readNamePattern = args$readNamePattern)
  }
  
  if(args$singleBarcode == FALSE | args$singleBarcode == "index2"){
    I2.parsed <- parLapply(
      cluster, 
      unique(samples.df$barcode2), 
      parseIndexReads,
      indexFilePath = args$index2,
      barcodeLength = args$barcode2Length,
      maxMismatch = args$maxMismatch,
      submat = submat,
      readNamePattern = args$readNamePattern)
  }
  
  stopCluster(cluster)
}else{
  if(args$singleBarcode == FALSE | args$singleBarcode == "index1"){
    I1.parsed <-  lapply(
      unique(samples.df$barcode1), 
      parseIndexReads,
      indexFilePath = args$index1,
      barcodeLength = args$barcode1Length,
      maxMismatch = args$maxMismatch,
      submat = submat,
      readNamePattern = args$readNamePattern)
  }
  
  if(args$singleBarcode == FALSE | args$singleBarcode == "index2"){
    I2.parsed <- lapply(
      unique(samples.df$barcode2), 
      parseIndexReads,
      indexFilePath = args$index2,
      barcodeLength = args$barcode2Length,
      maxMismatch = args$maxMismatch,
      submat = submat,
      readNamePattern = args$readNamePattern)
  }
}

if(args$singleBarcode == FALSE | args$singleBarcode == "index1"){
  names(I1.parsed) <- unique(samples.df$barcode1)
  pandoc.title("Barcode1 breakdown:")
  pandoc.table(sapply(I1.parsed, length), split.tables = Inf)
}

if(args$singleBarcode == FALSE | args$singleBarcode == "index2"){
  names(I2.parsed) <- unique(samples.df$barcode2)
  pandoc.title("Barcode2 breakdown:")
  pandoc.table(sapply(I2.parsed, length), split.tables = Inf)
}

if(args$singleBarcode == FALSE){
  demultiplexedReadNames <- mapply(
    function(barcode1, barcode2){
      intersect(I1.parsed[[barcode1]], I2.parsed[[barcode2]])
    },
    barcode1 = samples.df$barcode1,
    barcode2 = samples.df$barcode2)
  names(demultiplexedReadNames) <- paste0(
    samples.df$barcode1, samples.df$barcode2)
}else{
  if(args$singleBarcode == "index1"){
    demultiplexedReadNames <- I1.parsed
  }else{
    demultiplexedReadNames <- I2.parsed
  }
}

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
pandoc.title("Read counts for each sample.")
pandoc.table(samples.df, split.tables = Inf)
# Ambiguous reads
message(paste0("Ambiguous reads: ", length(ambiguousReads)))
# Unassigned reads
readNames <- id(readFastq(args$read1))
readNames <- str_extract(as.character(readNames), args$readNamePattern)
unassignedReadNames <- readNames[
  !readNames %in% unlist(demultiplexedReadNames)]
unassignedReadNames <- unassignedReadNames[
  !unassignedReadNames %in% ambiguousReads]
message(paste0("\nUnassigned reads: ", length(unassignedReadNames)))

# Create multiplex dataframe for subseting sequencing files
multiplexedData <- data.frame(
  "sampleName" = Rle(
    values = samples.df$sampleName, 
    length = sapply(demultiplexedReadNames, length)),
  "readName" = unlist(demultiplexedReadNames),
  row.names = NULL)

ambiguousData <- data.frame(
  "sampleName" = rep("ambiguous", length(ambiguousReads)),
  "readName" = ambiguousReads,
  row.names = NULL)

unassignedData <- data.frame(
  "sampleName" = rep("unassigned", length(unassignedReadNames)),
  "readName" = unassignedReadNames,
  row.names = NULL)

multiplexedData <- rbind(multiplexedData, ambiguousData, unassignedData)

stopifnot(nrow(multiplexedData) == length(readNames))

if(args$poolreps){
  multiplexedData$sampleName <- gsub("-\\d+$", "", multiplexedData$sampleName)
}

message(paste0("\nReads to be written to files: ", nrow(multiplexedData)))

# Write files to read files to outfolder directory
writeDemultiplexedSequences <- function(readFilePath, type, multiplexedData, 
                                        readNamePattern, outfolder, compress){
  # Require packages for parallel processing
  dependencies <- c("stringr", "ShortRead", "Biostrings")
  loaded <- sapply(dependencies, require, character.only = TRUE)
  stopifnot(all(loaded))
  # Load read sequences and sequence names then write to file
  reads <- readFastq(readFilePath)
  ids <- str_extract(as.character(id(reads)), readNamePattern)
  reads <- reads[match(multiplexedData$readName, ids)]
  reads <- split(reads, multiplexedData$sampleName)
  null <- lapply(1:length(reads), function(i, reads, type, outfolder, compress){
    if(compress){  
      filePath <- file.path(
        outfolder, paste0(names(reads[i]), ".", type, ".fastq.gz"))
    }else{
      filePath <- file.path(
        outfolder, paste0(names(reads[i]), ".", type, ".fastq"))
    }
    writeFastq(reads[[i]], file = filePath, compress = compress)
    message(
      paste0("\nWrote ", length(reads[[i]]), " reads to:\n", filePath, "."))
  }, reads = reads, type = type, outfolder = outfolder, compress = compress)
  return(list(readFilePath, type, outfolder))
}

if(args$cores > 0){
  cluster <- makeCluster(min(c(detectCores(), args$cores)))
  
  if(args$singleBarcode == FALSE){
    readList <- list("i1", "i2", "r1", "r2")
  }else if(args$singleBarcode == "index1"){
    readList <- list("i1", "r1", "r2")
  }else{
    readList <- list("i2", "r1", "r2")
  }
  
  if(args$singleBarcode == FALSE){
    readPaths <- list(args$index1, args$index2, args$read1, args$read2)
  }else if(args$singleBarcode == "index1"){
    readPaths <- list(args$index1, args$read1, args$read2)
  }else{
    readPaths<- list(args$index2, args$read1, args$read2)
  }
  
  demultiplex <- clusterMap(
    cluster,
    writeDemultiplexedSequences,
    readFilePath = readPaths,
    type = readList,
    MoreArgs = list(
      multiplexedData = multiplexedData,
      readNamePattern = args$readNamePattern,
      outfolder = args$outfolder,
      compress = args$compress))
  
  stopCluster(cluster)
}else{
  if(args$singleBarcode == FALSE){
    readList <- list("i1", "i2", "r1", "r2")
  }else if(args$singleBarcode == "index1"){
    readList <- list("i1", "r1", "r2")
  }else{
    readList <- list("i2", "r1", "r2")
  }
  
  if(args$singleBarcode == FALSE){
    readPaths <- list(args$index1, args$index2, args$read1, args$read2)
  }else if(args$singleBarcode == "index1"){
    readPaths <- list(args$index1, args$read1, args$read2)
  }else{
    readPaths<- list(args$index2, args$read1, args$read2)
  }
  
  demultiplex <- mapply(
    writeDemultiplexedSequences,
    readFilePath = readPaths,
    type = readList,
    MoreArgs = list(
      multiplexedData = multiplexedData,
      readNamePattern = args$readNamePattern,
      outfolder = args$outfolder,
      compress = args$compress))
}

message("\nDemultiplexing complete.")
q()
