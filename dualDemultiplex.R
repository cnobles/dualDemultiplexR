#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions
# Set Global options and load intiial packages ---------------------------------
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)

code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))

# Set up and gather command line arguments -------------------------------------
## Argument parser =============================================================
parser <- ArgumentParser(
  description = "R-based demultiplexing for dual-barcoded Illumina sequencing runs.")
parser$add_argument(
  "-m", "--manifest", type = "character", help = "Manifest file (*.yaml, *.yml, *.csv, *.tsv).")
parser$add_argument(
  "--read1", type = "character", default = "NA", help = "Path to Illumina R1 file (FASTQ).")
parser$add_argument(
  "--read2", type = "character", default = "NA", help = "Path to Illumina R2 file (FASTQ).")
parser$add_argument(
  "--index1", type = "character", default = "NA", help = "Path to Illumina I1 file (FASTQ).")
parser$add_argument(
  "--index2", type = "character", default = "NA", help = "Path to Illumina I2 file (FASTQ).")
parser$add_argument(
  "-o", "--outfolder", nargs = 1, type = "character", help = "Output folder.")
parser$add_argument(
  "-p", "--poolreps", action = "store_true", help = "Pools replicates.")
parser$add_argument(
  "--singleBarcode", action = "store_true",
  help = "Demultiplex with only a single barcode.")
parser$add_argument(
  "--barcode1", nargs = 1, type = "character", default = "I1", 
  help = "Read type containing barcode 1 sequences. (Default = 'I1', options include 'R1', 'R2', 'I1', and 'I2')")
parser$add_argument(
  "--barcode2", nargs = 1, type = "character", default = "I2", 
  help = "Read type containing barcode 2 sequences. (Default is 'I2', same options as --barcode1)")
parser$add_argument(
  "--barcode1Length", nargs = 1, type = "integer", default = 8, 
  help = "Length of barcode1, in nucleotides. Default = 8.")
parser$add_argument(
  "--barcode2Length", nargs = 1, type = "integer", default = 8,
  help = "Length of barcode2, in nucleotides. Default = 8.")
parser$add_argument(
  "--maxMismatch", nargs = 1, type = "integer", 
  help = "Max mismatch allowed in barcodes (Default = 0). It is recommended to use either ambiguous nucleotide codes or designate mismatch allowance, not both.")
parser$add_argument(
  "--bc1Mismatch", nargs = 1, type = "integer", default = 0,
  help = "Max mismatch allowed for barcode1 (Default = 0).")
parser$add_argument(
  "--bc2Mismatch", nargs = 1, type = "integer", default = 0,
  help = "Max mismatch allowed for barcode2 (Default = 0).")
parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", default = "[\\w:-]+",
  help = "Regular expression for pattern matching read names. Should not contain R1/R2/I1/I2 specific components. Default is [\\w:-]+")
parser$add_argument(
  "--compress", action = "store_true", help = "Output fastq files are gzipped.")
parser$add_argument(
  "-c", "--cores", nargs = 1, default = 0, type = "integer", 
  help = "Max cores to be used. If 0 (default), program will not utilize parallel processing.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

demulti <- data.frame(
  "readType" = c("R1", "R2", "I1", "I2"),
  "path" = c(args$read1, args$read2, args$index1, args$index2))
demulti$barcode1 <- grepl(args$barcode1, demulti$readType)
demulti$barcode2 <- grepl(args$barcode2, demulti$readType)

if(demulti$readType[demulti$barcode1] == demulti$readType[demulti$barcode2]){
  stop("Please select different read types for barcodes 1 and 2.")}

if(demulti$readType[demulti$barcode1] == "NA"){
  stop("Barcode 1 is set to a read type that is not provided.")}

if(demulti$readType[demulti$barcode2] == "NA"){
  stop("Barcode 2 is set to a read type that is not provided.")}

if(args$singleBarcode){
  demulti$barcode2 <- FALSE
}

if(!is.null(args$maxMisMatch)){
  args$bc1Mismatch <- args$maxMismatch
  args$bc2Mismatch <- args$maxMismatch
}

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("manifest :", "index1 :", "index2 :", "read1 :", "read2 :", 
          "outfolder :", "poolreps :", "singleBarcode :", "cores :", 
          "barcode1 :", "barcode2 :", "barcode1Length :", "barcode2Length :", 
          "bc1Mismatch :", "bc2Mismatch :", "readNamePattern :"),
        input_table$Variables),]
pandoc.title("Demultiplex Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Create output directory if not currently available ---------------------------
if(!file.exists(args$outfolder)){
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.")
}

# Load additional packages -----------------------------------------------------
add_packs <- c("stringr", "ShortRead", "Biostrings")
add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE))
if(!all(add_packs_loaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(add_packs_loaded), 
    "Loaded" = add_packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

source(file.path(
  code_dir, "supporting_scripts/binary_ambiguous_nucleotide_scoring_matrix.R"))

submat <- banmat()

# Load manifest / sample mapping file ------------------------------------------
fileExt <- unlist(strsplit(args$manifest, "\\."))
fileExt <- fileExt[length(fileExt)]

if(fileExt %in% c("yaml", "yml")){
  suppressMessages(library("yaml"))
  if(!"package:yaml" %in% search()) stop("Package:yaml not loaded or installed.")
  manifest <- yaml.load_file(args$manifest)
  if(args$singleBarcode){
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "barcode1" = sapply(manifest$samples, function(x) x$barcode1),
      row.names = NULL)
  }else{
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "barcode1" = sapply(manifest$samples, function(x) x$barcode1),
      "barcode2" = sapply(manifest$samples, function(x) x$barcode2),
      row.names = NULL)
  }
}else{
  if(fileExt == "csv"){
    manifest <- read.csv(args$manifest)
  }else if(fileExt == "tsv"){
    manifest <- read.delim(args$manifest)
  }
  if(args$singleBarcode){
    samples_df <- manifest[, c("sampleName", "barcode1")]
  }else{
    samples_df <- manifest[, c("sampleName", "barcode1", "barcode2")]
  }
}

if(!args$singleBarcode){
  uniqueSamples <- nrow(samples_df[,c("barcode1", "barcode2")]) == 
    nrow(unique(samples_df[,c("barcode1", "barcode2")]))
  if(!uniqueSamples) stop("Ambiguous barcoding of samples. Please correct.")
}else{
  uniqueSamples <- length(samples_df[,c("barcode1")]) == 
    length(unique(samples_df[,"barcode1"]))
  if(!uniqueSamples) stop("Ambiguous barcoding of samples. Please correct.")
}

# Read in barcode sequences ----------------------------------------------------
parseIndexReads <- function(barcode, indexFilePath, barcodeLength, maxMismatch, 
                            submat, readNamePattern){
  # Require packages for parallel processing
  dependencies <- c("stringr", "ShortRead", "Biostrings")
  loaded <- sapply(dependencies, require, character.only = TRUE)
  stopifnot(all(loaded))

  # Load index file sequences and sequence names
  index <- readFastq(indexFilePath)
  indexReads <- DNAStringSet(sread(index), start = 1, end = barcodeLength)
  names(indexReads) <- str_extract(as.character(id(index)), readNamePattern)
  
  # Trim barcode if necessary
  barcode <- as.character(DNAStringSet(barcode, start = 1, end = barcodeLength))

  # Identify read names with sequences above or equal to the minscore
  names(unlist(vmatchPattern(
    barcode, indexReads, max.mismatch = maxMismatch, fixed = FALSE)))
}

if(args$cores > 0){
  suppressMessages(library(parallel))
  cluster <- makeCluster(min(c(detectCores(), args$cores)))
  
  BC1_parsed <-  parLapply(
    cluster,
    unique(samples_df$barcode1), 
    parseIndexReads,
    indexFilePath = demulti$path[demulti$barcode1],
    barcodeLength = args$barcode1Length,
    maxMismatch = args$bc1Mismatch,
    submat = submat,
    readNamePattern = args$readNamePattern)
  
  names(BC1_parsed) <- unique(samples_df$barcode1)
  pandoc.title("Barcode1 breakdown:")
  pandoc.table(data.frame(
      "Barcode1" = names(BC1_parsed),
      "Read Counts" = sapply(BC1_parsed, length),
      row.names = NULL))
  
  if(!args$singleBarcode){
    BC2_parsed <- parLapply(
      cluster, 
      unique(samples_df$barcode2), 
      parseIndexReads,
      indexFilePath = demulti$path[demulti$barcode2],
      barcodeLength = args$barcode2Length,
      maxMismatch = args$bc2Mismatch,
      submat = submat,
      readNamePattern = args$readNamePattern)
  }
  
  stopCluster(cluster)
}else{
  BC1_parsed <-  lapply(
    unique(samples_df$barcode1), 
    parseIndexReads,
    indexFilePath = demulti$path[demulti$barcode1],
    barcodeLength = args$barcode1Length,
    maxMismatch = args$bc1Mismatch,
    submat = submat,
    readNamePattern = args$readNamePattern)

  names(BC1_parsed) <- unique(samples_df$barcode1)
  pandoc.title("Barcode1 breakdown:")
  pandoc.table(data.frame(
    "Barcode1" = names(BC1_parsed),
    "Read Counts" = sapply(BC1_parsed, length),
    row.names = NULL))
  
  if(!args$singleBarcode){
    BC2_parsed <- lapply(
      unique(samples_df$barcode2), 
      parseIndexReads,
      indexFilePath = demulti$path[demulti$barcode2],
      barcodeLength = args$barcode2Length,
      maxMismatch = args$bc2Mismatch,
      submat = submat,
      readNamePattern = args$readNamePattern)
  }
}

if(!args$singleBarcode){
  names(BC2_parsed) <- unique(samples_df$barcode2)
  pandoc.title("Barcode2 breakdown:")
  pandoc.table(data.frame(
    "Barcode2" = names(BC2_parsed),
    "Read Counts" = sapply(BC2_parsed, length),
    row.names = NULL))
}

if(!args$singleBarcode){
  demultiplexedReadNames <- mapply(
    function(barcode1, barcode2){
      intersect(BC1_parsed[[barcode1]], BC2_parsed[[barcode2]])
    },
    barcode1 = samples_df$barcode1,
    barcode2 = samples_df$barcode2,
    SIMPLIFY = FALSE)
  names(demultiplexedReadNames) <- paste0(
    samples_df$barcode1, samples_df$barcode2)
}else{
  demultiplexedReadNames <- BC1_parsed
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
samples_df$read_counts <- sapply(demultiplexedReadNames, length)
pandoc.title("Read counts for each sample.")
pandoc.table(samples_df, split.tables = Inf)
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

# Create multiplex dataframe for subseting sequencing files --------------------
multiplexedData <- data.frame(
  "sampleName" = Rle(
    values = samples_df$sampleName, 
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
    if(file.exists(filePath)) unlink(filePath)
    writeFastq(reads[[i]], file = filePath, compress = compress)
    message(
      paste0("\nWrote ", length(reads[[i]]), " reads to:\n", filePath, "."))
  }, reads = reads, type = type, outfolder = outfolder, compress = compress)
  return(list(readFilePath, type, outfolder))
}

if(args$cores > 0){
  cluster <- makeCluster(min(c(detectCores(), args$cores)))
  
  readList <- demulti$readType[demulti$path != "NA"]
  readPaths <- demulti$path[match(readList, demulti$readType)]

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
  readList <- demulti$readType[demulti$path != "NA"]
  readPaths <- demulti$path[match(readList, demulti$readType)]
  
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
