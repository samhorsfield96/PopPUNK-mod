import argparse
from multiprocessing import Pool
from KDE_distance import KDE_JS_divergence
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import itertools

def get_options():
    description = 'Comput pairwise 2D distances between selection of files'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python pairwise_2D_distances.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    required=True,
                    help='Directory containing run files from gridsearch run.')
    IO.add_argument('--outpref',
                    default = "output",
                    help='Output prefix. Default = "output"')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1.')

    return parser.parse_args()

def read_file(file):
    df = pd.read_csv(file, header=None, sep="\t", index_col=False)
    
    # downsample rows
    if len(df.columns) > 2:
        df = df[df.columns[[-2, -1]]]

    return df.to_numpy()

def get_distance(file_tuple):
    df1 = read_file(file_tuple[0])
    df2 = read_file(file_tuple[1])
    js_distance = KDE_JS_divergence(df1, df2, gamma=0.25, eps=0.0, log=False)
    min_core = np.min(df2[:, 0])
    max_core = np.max(df2[:, 0])
    min_acc = np.min(df2[:, 1])
    max_acc = np.max(df2[:, 1])
    return (file_tuple[0], file_tuple[1], js_distance, (min_core, max_core, min_acc, max_acc))

def main():
    options = get_options()

    indir = options.indir
    outpref = options.outpref
    threads = options.threads

    files = []
    for path in Path(indir).glob("*.txt"):
        # Print the path (file or directory) to the console
        files.append(str(path))
    
    #print(files)
    file_combinations = list(itertools.combinations(files, 2))

    with open(outpref + ".tsv", "w") as o:
        o.write("File1\tFile2\tDistance\n")
        with Pool(processes=threads) as pool:
            for file1, file2, js_distance, bounds in tqdm(pool.imap(get_distance, file_combinations), total = len(file_combinations)):
                o.write(f"{file1}\t{file2}\t{str(js_distance)}\n")

if __name__ == "__main__":
    main()