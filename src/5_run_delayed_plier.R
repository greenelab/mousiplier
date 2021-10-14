library(HDF5Array)
library(vroom)
library(dplyr)

source('delayed_plier.R')
setAutoRealizationBackend("HDF5Array") #supportedRealizationBackends(), getRealizationBackend()

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

# Load pathways -----------------------------------------------------------------------------------
pathways <- vroom('../data/plier_pathways.tsv', delim='\t')
pathway_genes <- pathways['...1']
colnames(pathway_genes) <- 'genes'
pathway_numbers <- subset(pathways, select=-c(1))
pathway_matrix <- as.matrix(pathway_numbers)
rownames(pathway_matrix) <- pathway_genes$genes

# Load counts -------------------------------------------------------------------------------------
expression_array <- DelayedArray(seed=HDF5ArraySeed(filepath="../data/no_scrna_tpm.h5", name="counts"))

# Get the row and column names  -------------------------------------------------------------------
expression_file <- '../data/no_scrna_tpm.tsv'

# This tibble uses about 12G of memory
expression_df <- vroom::vroom(expression_file, delim='\t', )

genes <- colnames(expression_df)
samples <- expression_df$sample
rm(expression_df)

# Set the row and column names for the expression array -------------------------------------------
# TODO add last sample to array and remove the indexing on this line
rownames(expression_array) <- samples[-length(samples)]
# Dear future me, negative indexing in R is different from negative indexing in Python. This says
# "get all items except the first one", not "get the last item"
colnames(expression_array) <- genes[-1]  # First entry is 'sample', so remove it

# Run PLIER ---------------------------------------------------------------------------------------
ptm <- proc.time()
PLIER.res <- PLIER(expression_array, pathway_matrix, output_path = "../data/")
print(proc.time()-ptm)
