"""
PLIER uses singular vectors as a starting point for its optimization. This script uses
incremental PCA to calculate PCs to use as a starting point without running out of memory
"""

import argparse
import os

import numpy as np
import pandas as pd
from sklearn.decomposition import IncrementalPCA
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('expression_file',
                        help='The tsv formatted expression file '
                             'produced by 3_preprocess_expression.py')
    parser.add_argument('out_dir', help='The directory to save the results to')
    parser.add_argument('--n_components',
                        help='The number of components to return from PCA',
                        default=200)
    args = parser.parse_args()

    # TODO test with different sizes to see if algo has converged

    pca = IncrementalPCA(n_components=args.n_components)

    columns_to_skip = 'sample'

    filter_fn = lambda x: x not in columns_to_skip

    CHUNKSIZE = 1000

    with pd.read_csv(args.expression_file,
                     chunksize=CHUNKSIZE,
                     delimiter='\t',
                     usecols=filter_fn) as reader:
        for chunk in tqdm(reader, total=190000 // CHUNKSIZE):
            data = chunk.to_numpy()
            if len(chunk) < args.n_components:
                continue
            pca.partial_fit(data)


    d = pca.singular_values_
    U = pca.components_.T

    transformed_chunks = []

    with pd.read_csv(args.expression_file,
                     chunksize=CHUNKSIZE,
                     delimiter='\t',
                     usecols=filter_fn) as reader:
        for i, chunk in tqdm(enumerate(reader), total=190000 // CHUNKSIZE):
            arr = chunk.to_numpy()

            # [samples x genes] x [genes x LVs] = [samples x LVs]
            transformed_chunk = arr @ U
            transformed_chunks.append(transformed_chunk)

    V = np.concatenate(transformed_chunks).T

    # Store results
    np.savetxt(os.path.join(args.out_dir, 'd.tsv'), d, delimiter='\t')
    np.savetxt(os.path.join(args.out_dir, 'U.tsv'), U, delimiter='\t')
    np.savetxt(os.path.join(args.out_dir, 'V.tsv'), V, delimiter='\t')
