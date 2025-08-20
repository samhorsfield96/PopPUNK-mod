import numpy as np
import pandas as pd
import argparse

def get_options():
    description = 'Return JC adjusted core genome.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python get_core_mu.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='Input distance file.')
    
    return parser.parse_args()

def read_distfile(filename):
    # read first line, determine if csv
    with open(filename, "r") as f:
        first_line = f.readline()
        if "," in first_line:
            obs = pd.read_csv(filename, index_col=None, header=None, sep=",")
        else:
            obs = pd.read_csv(filename, index_col=None, header=None, sep="\t")

    if len(obs.columns) == 2:
        obs.rename(columns={obs.columns[0]: "Core",
                           obs.columns[1]: "Accessory"}, inplace=True)
    elif len(obs.columns) == 4:
        # rename columns
        obs.rename(columns={obs.columns[0]: "Sample1", obs.columns[1] : "Sample2", obs.columns[2]: "Core",
                           obs.columns[3]: "Accessory"}, inplace=True)
    else:
        print("Incorrect number of columns in distfile. Should be 2 or 4.")
        sys.exit(1)

    obs['Core'] = pd.to_numeric(obs['Core'])
    obs['Accessory'] = pd.to_numeric(obs['Accessory'])

    obs = obs[['Core','Accessory']].to_numpy()

    return obs

def main():
    options = get_options()
    infile = options.infile

    obs_df = read_distfile(infile)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(np.max(obs_df[:,0]))
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))
    core_mu = max_real_core
    print("core_mu set to: {}".format(core_mu))

if __name__ == "__main__":
    main()