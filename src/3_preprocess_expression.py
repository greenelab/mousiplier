"""
This script converts counts to RPKM, row normalizes, and maps gene symbols for
a recount compendium
"""
from typing import Dict, Set

import argparse
import numpy as np
import tqdm


from utils import get_ensembl_mappings


def parse_gene_lengths(file_path: str) -> Dict[str, int]:
    """Parses a tsv file containing genes and their length

    Arguments
    ---------
    file_path - The path to the file mapping genes to lengths

    Returns
    -------
    gene_to_len - A dict mapping ensembl gene ids to their length in base pairs
    """
    gene_to_len = {}
    with open(file_path) as in_file:
        # Throw out header
        in_file.readline()
        for line in in_file:
            line = line.replace('"', '')
            gene, length = line.strip().split('\t')
            try:
                gene_to_len[gene] = int(length)
            except ValueError:
                # Some genes have no length, but will be removed in a later step
                pass
    return gene_to_len


def get_pathway_genes(pathway_file: str) -> Set[str]:
    """
    Read which genes are present in the pathway matrix file

    Arguments
    ---------
    pathway_file: The path to the file storing the pathway matrix as a genes x pathways tsv

    Returns
    -------
    pathway_genes: The set of all genes used in pathways
    """
    with open(args.pathway_file) as pathway_file:
        pathway_genes = set()

        header = pathway_file.readline()
        for line in pathway_file:
            line = line.strip().split('\t')
            gene = line[0]
            pathway_genes.add(gene)
        return pathway_genes


def calculate_tpm(counts: np.ndarray, gene_length_arr: np.ndarray) -> np.ndarray:
    """"Given an array of counts, calculate the transcripts per kilobase million
    based on the steps here:
    https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

    Arguments
    ---------
    counts: The array of transcript counts per gene
    gene_length_arr: The array of lengths for each gene in counts

    Returns
    -------
    tpm: The tpm normalized expression data
    """
    counts = np.array(counts, dtype=float)

    reads_per_kb = counts / gene_length_arr

    sample_total_counts = np.sum(counts)
    per_million_transcripts = sample_total_counts / 1e6

    tpm = reads_per_kb / per_million_transcripts

    return tpm


LINES_IN_FILE = 190112

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('count_file', help='The file containing the count matrix generated by '
                                           'remove_scrnaseq.py')
    parser.add_argument('gene_file', help='The file with gene lengths from get_gene_lengths.R')
    parser.add_argument('pathway_file', help='The file mapping genes to pathways')
    parser.add_argument('out_file', help='The file to save the normalized results to')

    args = parser.parse_args()

    # Map Ensembl to genesymbol
    ensembl_to_genesymbol = get_ensembl_mappings()

    # Get gene lengths to allow RPKM normalization
    gene_to_len = parse_gene_lengths(args.gene_file)

    pathway_genes = get_pathway_genes(args.pathway_file)

    # TPM normalize data
    with open(args.count_file, 'r') as count_file:
        header = count_file.readline()
        header = header.replace('"', '')
        header_genes = header.strip().split('\t')
        header_genes = [gene.split('.')[0] for gene in header_genes]

        header_gene_symbols = []
        for gene in header_genes:
            if gene in ensembl_to_genesymbol:
                header_gene_symbols.append(ensembl_to_genesymbol[gene])
            else:
                header_gene_symbols.append(None)

        bad_indices = []
        # Keep only the first instance of each gene in the case that multiple
        # Ensembl genes get mapped to one gene symbol

        duplicate_count = 0
        not_pathway_count = 0
        no_length_count = 0
        genes_seen = set()
        for i, gene in enumerate(header_gene_symbols):
            if gene is None or gene in genes_seen:
                bad_indices.append(i)
            # Remove genes that aren't in our prior pathways
            elif gene not in pathway_genes:
                bad_indices.append(i)
            else:
                genes_seen.add(gene)

        # Remove genes with unknown lengths
        gene_length_arr = []
        for i, gene in enumerate(header_genes):

            if gene not in gene_to_len.keys():
                bad_indices.append(i)
                gene_length_arr.append(None)
            else:
                gene_length_arr.append(gene_to_len[gene])

        # sort bad_indices and deduplicate
        bad_indices = list(set(bad_indices))
        bad_indices.sort()

        for index in reversed(bad_indices):
            del gene_length_arr[index]
        gene_length_arr = np.array(gene_length_arr)

        means = None
        M2 = None
        maximums = None
        minimums = None

        samples_seen = set()
        # First time through the data, calculate statistics
        for i, line in tqdm.tqdm(enumerate(count_file), total=LINES_IN_FILE):
            line = line.replace('"', '')
            line = line.strip().split('\t')
            sample = line[0]

            # Remove duplicates
            if sample in samples_seen:
                continue
            samples_seen.add(sample)

            try:
                # Thanks to stackoverflow for this smart optimization
                # https://stackoverflow.com/a/11303234/10930590

                counts = line[1:]  # bad_indices is still correct because of how R saves tables
                for index in reversed(bad_indices):
                    del counts[index]

                tpm = calculate_tpm(counts, gene_length_arr)

                if any(np.isnan(tpm)):
                    continue

                # Online variance calculation https://stackoverflow.com/a/15638726/10930590
                if means is None:
                    means = tpm
                    M2 = 0
                    maximums = tpm
                    minimums = tpm
                else:
                    delta = tpm - means
                    means = means + delta / (i + 1)
                    M2 = M2 + delta * (tpm - means)
                    maximums = np.maximum(maximums, tpm)
                    minimums = np.minimum(minimums, tpm)

            except ValueError as e:
                # Throw out malformed lines caused by issues with downloading data
                print(e)

        per_gene_variances = M2 / (i-1)
        max_min_diff = maximums - minimums

        n = args.num_genes

        out_file = open(args.out_file, 'w')

        header = header.strip().split('\t')

        # Use numpy to allow indexing with a list of indices
        header_arr = np.array(header)
        index_mask = np.ones(header_arr.shape)
        index_mask[bad_indices] = False
        header = header_arr[index_mask]

        header = header_arr.tolist()

        header = 'sample\t' + '\t'.join(header)
        out_file.write(header)
        out_file.write('\n')

    with open(args.count_file, 'r') as count_file:
        samples_seen = set()
        # Second time through the data - standardize and write outputs
        for i, line in tqdm.tqdm(enumerate(count_file), total=LINES_IN_FILE):
            line = line.replace('"', '')
            line = line.strip().split('\t')
            sample = line[0]

            if sample in samples_seen:
                continue
            samples_seen.add(sample)

            try:
                counts = line[1:]  # bad_indices is still correct because of how R saves tables
                for index in reversed(bad_indices):
                    del counts[index]

                tpm = calculate_tpm(counts, gene_length_arr)

                if any(np.isnan(tpm)):
                    continue

                # Zero-one standardize
                standardized_tpm = (tpm - minimums) / max_min_diff

                # Keep only most variable genes
                tpm_list = standardized_tpm.tolist()
                tpm_strings = ['{}'.format(x) for x in tpm_list]

                out_file.write('{}\t'.format(sample))
                out_file.write('\t'.join(tpm_strings))
                out_file.write('\n')

            except ValueError as e:
                # Throw out malformed lines caused by issues with downloading data
                print(e)