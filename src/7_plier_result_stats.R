library(dplyr)
library(tibble)

# Make sure R's working directory is in the correct spot to make relative paths resolve correctly
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # If running as a script, finding the file is harder
  # https://stackoverflow.com/a/55322344/10930590
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)

  setwd(dirname(this_file))
}

plier_results <- readRDS('../output/plier.rds')

Z <- plier_results['Z'][[1]]
U <- plier_results['U'][[1]]
C <- plier_results['C'][[1]]
lambda <- plier_results['L2'][[1]]

# Write pathway histogram to file
png('../output/lv_per_pathway_hist.png', width=720, height=720)
# Using U != 0 because U isn't a binary matrix, it's just sparse due to L1 norm
hist(colSums(U != 0), main='Distribution of pathways per LV ', xlab='Pathways per latent variable')
dev.off()

# Calculate genes per pathway (C)
png('../output/gene_per_pathway_hist.png', width=720, height=720)
hist(colSums(C)[colSums(C) < 200], main='Distribution of genes per pathway (truncated at 200)',
     xlab='Genes per pathway', breaks=50)
dev.off()

png('../output/percent_genes_used_hist.png', width=720, height=720)
hist((colSums(Z != 0) / length(Z[,1])) * 100, breaks=20, main='Percent of all genes used per LV',
     xlab='Percent of genes used')
dev.off()

# Write results to tsvs
write.table(Z, file='../output/Z.tsv', sep='\t')
write.table(U, file='../output/U.tsv', sep='\t')
fileConn<-file('../output/lambda.txt')
writeLines(toString(lambda), fileConn)
close(fileConn)
