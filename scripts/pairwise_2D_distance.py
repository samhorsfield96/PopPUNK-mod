import argparse
from multiprocessing import Pool
from KDE_distance import KDE_JS_divergence
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import itertools
import sys
import math

def get_options():
    description = 'Comput pairwise 2D distances between selection of files'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python pairwise_2D_distances.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    default=None,
                    help='Directory containing run files from gridsearch run.')
    IO.add_argument('--infile',
                    default=None,
                    help='File generated from create_ppmod_input.py.')
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
    input1 = file_tuple[0]
    input2 = file_tuple[1]

    df1 = read_file(input1[0])
    df2 = read_file(input2[0])

    js_distance = KDE_JS_divergence(df1, df2, gamma=0.25, eps=0.0, log=False)
    min_core = np.min(df2[:, 0])
    max_core = np.max(df2[:, 0])
    min_acc = np.min(df2[:, 1])
    max_acc = np.max(df2[:, 1])

    if len(input1) == 1:
        distance = js_distance
    else:
        # calculate distance as done in pp-mod using euclidean distance
        js_distance1 = 0
        core1 = input1[1]
        inter1 = input1[2]
        rare1 = input1[3]
        results1 = (js_distance1, core1, inter1, rare1)

        js_distance2 = js_distance
        core2 = input2[1]
        inter2 = input2[2]
        rare2 = input2[3]
        results2 = (js_distance2, core2, inter2, rare2)

        euclidean_dist = math.dist(results1, results2)
        distance = euclidean_dist
    
    return (input1[0], input2[0], distance, (min_core, max_core, min_acc, max_acc))

def main():
    options = get_options()

    indir = options.indir
    infile = options.infile
    outpref = options.outpref
    threads = options.threads

    if indir == None and infile == None:
        print("Please specify one of --indir or --infile")
        sys.exit(1)

    files = []

    if indir != None:
        for path in Path(indir).glob("*.txt"):
            # Print the path (file or directory) to the console
            files.append((str(path)))
    else:
        with open(infile, "r") as f:
            f.readline()
            for line in f:
                split_line = line.rstrip().split("\t")
                pan_genes = int(split_line[1])
                core_genes = int(split_line[2]) / pan_genes
                intermediate_genes = int(split_line[3]) / pan_genes
                rare_genes = (split_line[4]) / pan_genes
                file = split_line[6]
                files.append((file, core_genes, intermediate_genes, rare_genes))

    #print(files)
    file_combinations = list(itertools.combinations(files, 2))

    with open(outpref + ".tsv", "w") as o:
        o.write("File1\tFile2\tDistance\n")
        with Pool(processes=threads) as pool:
            for file1, file2, js_distance, bounds in tqdm(pool.imap(get_distance, file_combinations), total = len(file_combinations)):
                o.write(f"{file1}\t{file2}\t{str(js_distance)}\n")

if __name__ == "__main__":
    main()