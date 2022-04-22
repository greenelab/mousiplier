"""
This script is used to reformat the count output from featureCount: i.e., combine
count from different samples
"""


import pandas as pd
import argparse

def get_name_map(infile):
    """
    :param infile: the map from bam file to sample ID
    example:   Cocaine_NAc_rep1Aligned.out.bam.sorted.bam	day1_cocaine_NAc_rep1
    :return: a dict containing the mapping information
    """
    name_map = {}
    with open(infile) as fh:
        for line in fh:
            bamName, id = line.strip().split()
            name_map[bamName] = id
    return name_map


def parse_arguments():
    """
    parse the command line argument
    :return: a parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("day1_count", help="the count input for day1 abstinence")
    parser.add_argument("day28_count", help="the count input for day28 abstinence")
    parser.add_argument("name_info", help="the map from bam file to sample ID")
    parser.add_argument("outfile", help="output file")

    return parser

def main():
    """
    :param argv: command line arguments
    the first is day1 count, and the second is day28 count, the third is name_map.txt
    the fourth is output file
    :return:
    """
    parser = parse_arguments()
    args = parser.parse_args()
    day1_counts = pd.read_csv(args.day1_count, header = 0, sep='\t', index_col = 0)
    day2_counts = pd.read_csv(args.day28_count, header = 0, sep='\t', index_col = 0)
    ### transpose so that rows are samples, and columns are genes
    day1_counts_tranposed = day1_counts.T
    day2_counts_tranposed = day2_counts.T
    name_map = get_name_map(args.name_info)
    outfh = open(args.outfile, 'w')

    if list(day1_counts_tranposed.columns.values) == list(day2_counts_tranposed.columns.values):
        # write header
        temp_line = ''
        for gene in list(day1_counts_tranposed.columns.values):
            temp_line += gene + '\t'
        temp_line = temp_line.strip() + '\n'
        outfh.write(temp_line)

        ### write for day1 samples
        for index, row in day1_counts_tranposed.iterrows():
            print("The index is: " + index)
            # write the sample
            outfh.write(name_map[index] + '\t')
            # write count
            temp_line = ''
            for c in row:
                temp_line += str(c) + '\t'
            temp_line = temp_line.strip() + '\n'
            outfh.write(temp_line)

        ### write for day28 samples
        for index, row in day2_counts_tranposed.iterrows():
            # write the sample
            outfh.write(name_map[index] + '\t')
            # write count
            temp_line = ''
            for c in row:
                temp_line += str(c) + '\t'
            temp_line = temp_line.strip() + '\n'
            outfh.write(temp_line)

    outfh.close()


if __name__ == "__main__":
    main()
