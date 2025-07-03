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
    IO.add_argument('--core-min-bounds',
                    type=str,
                    default = None,
                    help='Min value proportion bounds for core values for simulation to be chosen. Must be comma separated as 0.0,0.1')
    IO.add_argument('--core-max-bounds',
                    type=str,
                    default = None,
                    help='Max value proportion bounds for core values for simulation to be chosen. Must be comma separated as 0.0,0.1')
    IO.add_argument('--acc-min-bounds',
                    type=str,
                    default = None,
                    help='Min value proportion bounds for accessory values for simulation to be chosen. Must be comma separated as 0.0,0.1')
    IO.add_argument('--acc-max-bounds',
                    type=str,
                    default = None,
                    help='Max value proportion bounds for accessory values for simulation to be chosen. Must be comma separated as 0.0,0.1')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1.')            

    return parser.parse_args()

def get_distance(file, df1):
    df2 = np.loadtxt(file, delimiter='\t', dtype='float64')
    js_distance = KDE_JS_divergence(df1, df2, 1e-12, log=True)
    min_core = np.min(df2[:, 0])
    max_core = np.max(df2[:, 0])
    min_acc = np.min(df2[:, 1])
    max_acc = np.max(df2[:, 1])
    return (file, js_distance, (min_core, max_core, min_acc, max_acc))

def check_bounds(value, bounds):
    in_range = False
    if bounds[0] <= value <= bounds[1]:
        in_range = True

    return in_range

def main():
    options = get_options()
    indir = options.indir
    ref = options.ref
    outpref = options.outpref
    threads = options.threads
    core_min_bounds = options.core_min_bounds
    acc_min_bounds = options.acc_min_bounds
    core_max_bounds = options.core_max_bounds
    acc_max_bounds = options.acc_max_bounds

    files = []
    for path in Path(indir).glob("**/*.tsv"):
        # Print the path (file or directory) to the console
        files.append(str(path))
    
    ref_df = pd.read_csv(ref, header=None, sep="\t", index_col=False)
    
    # downsample rows
    if len(ref_df.columns) > 2:
        ref_df = ref_df[ref_df.columns[[-2, -1]]]

    ref_df = ref_df.to_numpy()
    #print(ref_df)

    max_core_ref = np.max(ref_df[:, 0])
    max_acc_ref = np.max(ref_df[:, 1])
    print(f"max_core_ref: {max_core_ref}, max_acc_ref: {max_acc_ref}")

    if core_min_bounds is not None:
        core_min_bounds = [float(x) * max_core_ref for x in core_min_bounds.split(",")]

    if core_max_bounds is not None:
        core_max_bounds = [float(x) * max_core_ref for x in core_max_bounds.split(",")]

    if acc_min_bounds is not None:
        acc_min_bounds = [float(x) * max_acc_ref for x in acc_min_bounds.split(",")]

    if acc_max_bounds is not None:
        acc_max_bounds = [float(x) * max_acc_ref for x in acc_max_bounds.split(",")]

    print(f"core_min_bounds: {core_min_bounds}, core_max_bounds: {core_max_bounds}")
    print(f"acc_min_bounds: {acc_min_bounds}, acc_max_bounds: {acc_max_bounds}")
    
    min_distance = (None, 1.0)
    with Pool(processes=threads) as pool:
        for file, js_distance, bounds in tqdm(pool.imap(partial(get_distance, df1=ref_df), files), total = len(files)):
            if js_distance < min_distance[1]:
                core_passed = True
                acc_passed = True
                min_core, max_core, min_acc, max_acc = bounds
                #print(bounds)
                # check bounds of selection
                if core_min_bounds is not None and core_max_bounds is not None:
                    if not (check_bounds(min_core, core_min_bounds) and check_bounds(max_core, core_max_bounds)):
                        core_passed = False
                if acc_min_bounds is not None and acc_max_bounds is not None:
                    if not (check_bounds(min_acc, acc_min_bounds) and check_bounds(max_acc, acc_max_bounds)):
                        core_passed = False

                if core_passed and acc_passed:
                    min_distance = (file, js_distance, bounds)
                # else:
                #     print(bounds)
                #     print(f"check_bounds(min_core, core_min_bounds): {check_bounds(min_core, core_min_bounds)}")
                #     print(f"check_bounds(max_core, core_max_bounds): {check_bounds(max_core, core_max_bounds)}")
                #     print(f"check_bounds(min_acc, acc_min_bounds): {check_bounds(min_acc, acc_min_bounds)}")
                #     print(f"check_bounds(max_acc, acc_max_bounds): {check_bounds(max_acc, acc_max_bounds)}")
    
    if len(min_distance) == 2:
        print("No best fitting distribution found.")
    else:
        print(f"File: {min_distance[0]}, Min distance: {min_distance[1]}, Min core: {min_distance[2][0]}, Max core: {min_distance[2][1]}, Min acc: {min_distance[2][2]}, Max acc: {min_distance[2][3]}")


if __name__ == "__main__":
    main()