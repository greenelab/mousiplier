import argparse


def parse_sample_files(file_paths):
    """
    Parse the list of samples to remove from the compendium
    """
    samples = []
    for file_path in file_paths:
        with open(file_path) as in_file:
            # Toss header
            in_file.readline()
            for line in in_file:
                sample = line.split(',')[0]
                samples.append(sample)
    return set(samples)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('compendium_counts', help='Recount3 compendium in count format')
    parser.add_argument('out_file', help='Path to save the results to')
    parser.add_argument('sample_files',
                        help='Metadata files from NCBI sample selector containing '
                             'samples to remove from compendium',
                        nargs='+')
    args = parser.parse_args()

    holdout_samples = parse_sample_files(args.sample_files)

    out_file = open(args.out_file, 'w')
    with open(args.compendium_counts) as in_file:
        header = in_file.readline()
        out_file.write(header)
        for line in in_file:
            sample = line.split('\t')[0]
            sample = sample.strip('"')
            if sample not in holdout_samples:
                out_file.write(line)
    out_file.close()