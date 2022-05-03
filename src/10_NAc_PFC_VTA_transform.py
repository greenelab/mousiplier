import pandas as pd
import argparse

from transform import PlierTransform

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('weight_file', help="weight info from Z loading")
    parser.add_argument('lambda_file', help="The file with lambda")
    parser.add_argument('expression_file', help="fpkm normalized expression data")
    parser.add_argument('outfile', help="The output file to save the values of latent vairable")
    args = parser.parse_args()

    ### read Z loading
    lv_df = pd.read_csv(args.weight_file, sep='\t')

    ### read expression data
    expression_df = pd.read_csv(args.expression_file, delimiter='\t', index_col=0)

    ### select the genes present in Z loading
    reformatted_expression_df = pd.DataFrame(index=expression_df.index)
    for i in lv_df.index:
        if i in expression_df.columns.values:
            reformatted_expression_df[i] = expression_df[i]
        else:
            reformatted_expression_df[i] = [0]*expression_df.shape[0]

    ### transform the gene expression into latent space
    transformer = PlierTransform(args.weight_file, args.lambda_file)
    transformed_df = transformer.transform(reformatted_expression_df)
    transformed_df.to_csv(args.outfile, sep="\t")
