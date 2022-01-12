library(vroom)
library(dplyr)
library(PLIER)

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

# Get the row and column names  -------------------------------------------------------------------
expression_file <- '../data/no_scrna_tpm.tsv'

# This tibble uses about 12G of memory
expression_df <- vroom::vroom(expression_file, delim='\t', )

genes <- colnames(expression_df)
samples <- expression_df$sample
expression_df <- expression_df[,-1]

expression_array <- as.matrix(expression_df)

rm(expression_df)

# Set the row and column names for the expression array -------------------------------------------
rownames(expression_array) <- samples
# Dear future me, negative indexing in R is different from negative indexing in Python. This says
# "get all items except the first one", not "get the last item"
colnames(expression_array) <- genes[-1]  # First entry is 'sample', so remove it

# Load PCA results for initializing PLIER ---------------------------------------------------------
d <- read.csv('../data/d.tsv', sep='\t', header=FALSE)
# Coerce d into a vector so `diff` works
d <- c(as.matrix(d))
V <- read.csv('../data/V.tsv', sep='\t', header=FALSE)

svdres <- list()
svdres$d <- d
# Transpose from 200 x 190k to 190k x 200
svdres$v <- t(V)

# Run PLIER ---------------------------------------------------------------------------------------
expression_array <- t(expression_array)
ptm <- proc.time()
PLIER.res <- PLIER(expression_array, pathway_matrix, minGenes=6, svdres=svdres)
saveRDS(PLIER.res, '../output/plier.rds')
print(proc.time()-ptm)
