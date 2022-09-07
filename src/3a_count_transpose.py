import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('original_count', help="original count file")
    parser.add_argument('transposed_count', help="transposed count file")
    args = parser.parse_args()

    counts = pd.read_csv(args.original_count, header=0, sep='\t', index_col=0)
    counts_T = counts.T
    counts_T.to_csv(args.transposed_count, sep='\t')
