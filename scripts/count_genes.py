import argparse
from collections import defaultdict
import statistics
import pandas as pd
import os

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Count genes and calculate statistics in terms of core, pangenome size, average gene frequency presence/absence matrices')

    # input options
    parser.add_argument('--infile', help='Path to presense/absence file',
                                    required=True)
    parser.add_argument('--header', help='Header present in input file',
                                    default=False,
                                    action="store_true")
    parser.add_argument('--core', help='Core gene frequency threshold (default = 0.95)',
                                    type=float,
                                    default=0.95)
    parser.add_argument('--rare', help='Rare gene frequency threshold (default = 0.05)',
                                    type=float,
                                    default=0.05)
    parser.add_argument('--min-freq', help='Minimum gene frequency to be counted. Can be integer (absolute count) or decimal (frequency). (default = 0.0)',
                                    type=float,
                                    default=0.0)

    return parser.parse_args()

def main():

    # Check input ok
    args = get_options()
    header = 0 if args.header == True else None

    basename = os.path.splitext(os.path.basename(args.infile))[0]

    if ".csv" in args.infile:
        df = pd.read_csv(args.infile, sep=",", header=header)
    elif ".tsv" in args.infile:
        df = pd.read_csv(args.infile, sep="\t", header=header)
    else:
        df = pd.read_csv(args.infile, sep=None, engine="python", header=header)

    n_genes = len(df.columns)
    n_indivuduals = len(df)

    # stats per individuals
    sum_rows = df.sum(axis=1)
    avg_genome_size = round(sum_rows.mean())
    std_genome_size = round(sum_rows.std())

    # core genome determination
    core_genome_threshold = n_indivuduals * args.core
    rare_genome_threshold = n_indivuduals * args.rare

    # stats per gene
    sum_cols = df.sum(axis=0)
    core_genome_size = (sum_cols >= core_genome_threshold).sum()
    rare_genome_size = ((sum_cols < rare_genome_threshold) & (sum_cols > 0)).sum()
    missing_genome_size = (sum_cols == 0).sum()
    intermediate_genome_size = ((sum_cols >= rare_genome_threshold) & (sum_cols < core_genome_threshold)).sum()

    print(f"{basename},core={core_genome_size},intermediate={intermediate_genome_size},rare={rare_genome_size},avg_genome_size={avg_genome_size},std_genome_size={std_genome_size}")

if __name__ == "__main__":
    main()
