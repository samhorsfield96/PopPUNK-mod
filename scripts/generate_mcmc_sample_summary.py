import argparse
import glob
from pathlib import Path
import pandas as pd
import os

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Generates summary of MCMC results per species')

    # input options
    parser.add_argument('--indir', help='Path to directory containing output from PopPUNK-mod outputs',
                                    required=True)
    parser.add_argument('--outpref', help='Output prefix',
                                    default="pp-mod_mcmc_summary.txt")
    parser.add_argument('--parse', default="",
                                    help='pattern to parse file names at to generate sample ID e.g. _core_mu_')

    return parser.parse_args()

# search recursively
def find_files_with_extension(root_dir, extension):
    root = Path(root_dir)
    yield from root.rglob(f'*{extension.lstrip(".")}')

def main():
    options = get_options()

    indir = options.indir
    outpref = options.outpref
    parse = options.parse

    rows = []
    count = 0
    for file in find_files_with_extension(indir, "_mcmc_posterior_samples.csv"):
        base = os.path.basename(file).split(parse)[0]
        df = pd.read_csv(file, header=0)
        median = df.quantile(0.5)
        lower_CI = df.quantile(0.025)
        upper_CI = df.quantile(0.975)
        
        # concatenate by column
        combined = pd.concat([median, lower_CI, upper_CI], axis=1)
        combined.columns = ["median", "lower_CI", "upper_CI"]

        # transpose and convert to wide format
        row = combined.stack().to_frame().T
        row.columns = [f"{param}_{stat}" for param, stat in row.columns]
        row.index = [base]
        row.index.name = "taxa"

        rows.append(row)

        count += 1
        if count % 10 == 0:
            print(f"At count: {count}")

    big_df = pd.concat(rows)
    big_df.index.name = "taxa"

    big_df.to_csv(outpref + ".csv")
        

if __name__ == "__main__":
    main()