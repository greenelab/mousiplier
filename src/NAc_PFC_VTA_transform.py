import pandas as pd
import sys

from transform import PlierTransform


(weight_file, lambda_file, expression_file, outfile) = sys.argv[1:]

### read Z loading
lv_df = pd.read_csv(weight_file, sep='\t')

### read expression data
expression_df = pd.read_csv(expression_file, delimiter='\t', index_col=0)

### select the genes present in Z loading
reformatted_expression_df = pd.DataFrame(index=expression_df.index)
for i in lv_df.index:
    if i in expression_df.columns.values:
        reformatted_expression_df[i] = expression_df[i]
    else:
        reformatted_expression_df[i] = [0]*expression_df.shape[0]
#reformatted_expression_df.to_csv(outfile, sep='\t')
#expression_df = pd.read_csv(sys.argv[1],delimiter='\t', index_col=0)

### transform the gene expression into latent space
transformer = PlierTransform(weight_file, lambda_file)
transformed_df = transformer.transform(reformtted_expression_df)
transformed_df.to_csv(outfile, sep="\t")
