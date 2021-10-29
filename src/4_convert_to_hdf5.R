# This script prepares recount3 data to be run through PLIER. It is based on
# Qiwen and Jackie's script here: https://github.com/greenelab/rheum-plier-data/blob/05822c3a7b2e7c663d79426b1117a36dc2a0b892/recount2/2-prep_recount_for_plier.R

library(tibble)
library(tidyr)
library(vroom)
library(dplyr)
library(biomaRt)
library(rhdf5)

read_header = function(file_path) {
  open_file = file(file_path, 'r')
  header = read.csv(open_file, nrows=1, header=TRUE, sep='\t')
  return(header)
}

append_h5 <- function(data, h5file, row_start, size) {
  row_extent <- seq_len(size) + row_start
  print(row_extent)
  ncols <- ncol(data)

  # Write the chunk of data to the correct location in the array
  h5write(data, h5file, 'counts', index=list(row_extent, seq_len(ncols)))
  h5closeAll()
}

processFile = function(filepath, h5_file, coltypes) {
  con = file(filepath, "r")
  # Throw away header
  readLines(con, n=1)
  i <- 0
  chunksize <- 5000
  while ( TRUE ) {
    chunk <- try(read.csv(con, nrows=chunksize, sep='\t', colClasses=coltypes, header=FALSE))
    # Drop sample column
    if ( length(chunk) == 0 ) {
      break
    }
    chunk[,1] <- NULL
    chunk <- as.matrix(chunk)

    append_h5(chunk, h5_file, i, NROW(chunk))
    
    i = i + NROW(chunk)
    
    # If we're out of data calling read.csv again will throw an error
    if (NROW(chunk) < chunksize){
      break
    }
  }
  close(con)
}

# Set directory -----------------------------------------------------------------------------------

# Make sure R's working directory is in the correct spot to make relative paths resolve correctly
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else{
  # If running as a script, finding the file is harder
  # https://stackoverflow.com/a/55322344/10930590
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  
  setwd(dirname(this_file))
}


# Load expression data ----------------------------------------------------------------------------
expression_file <- '../data/no_scrna_tpm.tsv'

# This tibble uses about 12G of memory
expression_df <- vroom::vroom(expression_file, delim='\t', )

genes <- colnames(expression_df)
samples <- expression_df$sample
rm(expression_df)

out_path <- '../data/no_scrna_tpm.h5'
if (!file.exists(out_path)){
  h5createFile(out_path)
  # length(genes-1) to ignore the 'sample' column in the header
  h5createDataset(out_path, 'counts', dims=c(length(samples), length(genes)-1))
}

coltypes = c('character', rep('numeric', times=length(genes)-1))

processFile(expression_file, out_path, coltypes)
