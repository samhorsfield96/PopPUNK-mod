import argparse
from multiprocessing import Pool
from functools import partial
from KDE_distance import KDE_JS_divergence
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

def get_options():
    description = 'Compare distributions from grid searchh to find optimal parameters'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python distance_to_gridsearch.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    required=True,
                    help='Directory containing run files from gridsearch run.')
    IO.add_argument('--ref',
                    required=True,
                    help='Reference file to be compared to.')
    IO.add_argument('--outpref',
                    default = "output",
                    help='Output prefix. Default = "output"')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1.')            

    return parser.parse_args()

def get_distance(file, df1):
    df2 = np.loadtxt(file, delimiter='\t', dtype='float64')
    js_distance = KDE_JS_divergence(df1, df2)
    return (file, js_distance)

def main():
    options = get_options()
    indir = options.indir
    ref = options.ref
    outpref = options.outpref
    threads = options.threads

    files = []
    for path in Path(indir).glob("**/*.tsv"):
        # Print the path (file or directory) to the console
        files.append(str(path))
    
    ref_df = pd.read_csv(ref, header=None, sep="\t", index_col=False)
    
    # downsample rows
    if len(ref_df.columns) > 2:
        ref_df = ref_df[ref_df.columns[[-2, -1]]]

    ref_df = ref_df.to_numpy()
    
    min_distance = (None, 1.0)
    with Pool(processes=threads) as pool:
        for file, js_distance in tqdm(pool.imap(partial(get_distance, df1=ref_df), files), total = len(files)):
            if js_distance < min_distance[1]:
                min_distance = (file, js_distance)
    
    print(f"File: {min_distance[0]}, Min distance: {min_distance[1]}")


if __name__ == "__main__":
    main()