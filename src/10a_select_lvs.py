import argparse

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('lv_file', help='File containing LVs to be filtered by variance')
    parser.add_argument('out_file', help='File to store the filtered LVs to ')
    parser.add_argument('--n_to_keep', help='The number of lvs to keep', default=10)
    args = parser.parse_args()

    lv_df = pd.read_csv(args.lv_file, sep='\t')
    samples = lv_df['sample']
    lv_df = lv_df.drop('sample', axis='columns')
    lv_df = lv_df.sample(frac=1).reset_index(drop=True)

    top_lvs = lv_df.var().nlargest(args.n_to_keep).index

    filtered_lvs = lv_df.loc[:,top_lvs]

    filtered_lvs['sample'] = samples
    filtered_lvs = filtered_lvs.set_index('sample')

    filtered_lvs.to_csv(args.out_file, sep='\t')