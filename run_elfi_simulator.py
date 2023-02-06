from run_elfi import gen_distances_elfi
from simulate_divergence import read_distfile
import numpy as np
import argparse
import os
import sys
from functools import partial
from multiprocessing import Pool
import tqdm
import matplotlib.pyplot as plt

def get_options():
    description = 'Fit model to PopPUNK data using Approximate Baysesian computation'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_ELFI.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--run-mode',
                    required=True,
                    choices=['sim', 'sample'],
                    help='Which run mode to specify. Choices are "sim" or "sample".')
    IO.add_argument('--core-size',
                    type=int,
                    default=10000,
                    help='Number of positions in core genome. Default = 10000 ')
    IO.add_argument('--pan-size',
                    type=int,
                    default=10000,
                    help='Number of positions in pangenome. Default = 10000 ')
    IO.add_argument('--pop-size',
                    type=int,
                    default=1000,
                    help='Population size for Wright-Fisher model. Default = 1000 ')
    IO.add_argument('--ngen',
                    type=int,
                    default=100,
                    help='Number of generations for Wright-Fisher model. Default = 100 ')
    IO.add_argument('--base-mu',
                    default="0.25,0.25,0.25,0.25",
                    help='Mutation rates from all other bases to each base, in order "A,C,G,T". Default = "0.25,0.25,0.25,0.25" ')
    IO.add_argument('--avg-gene-freq',
                    type=float,
                    default=0.5,
                    help='Average gene frequency in accessory genome. '
                         'Determines gene gain/loss rate e.g. 0.1 = gene gain/loss rate 1:9 '
                         'Default = "0.5" ')
    IO.add_argument('--distfile',
                    required=True,
                    help='popPUNK distance file to fit to. ')
    IO.add_argument('--params',
                    required=True,
                    help='Parameters file to run simulation.')
    IO.add_argument('--outpref',
                    default="sim",
                    help='Output prefix. Default = "sim"')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')

    return parser.parse_args()

def run_sim(index, params_list, max_hamming_core, max_jaccard_acc):
    param_set = params_list[index]
    size_pan = int(param_set[0])
    size_core = int(param_set[1])
    avg_gene_freq = float(param_set[2])
    base_mu = param_set[3]
    n_gen = int(param_set[4])
    pop_size = int(param_set[5])
    prop_gene = float(param_set[6])
    gene_gl = float(param_set[7])

    base_mu = [float(i) for i in base_mu.split(",")]

    # set evenly spaced core hamming values across generations
    core_mu = max_real_core / n_gen

    # round to 6 dp
    base_mu = [round(i, 6) for i in base_mu]

    # ensure probabilities sum to 1
    if sum(base_mu) != 1:
        base_mu[-1] = 1 - sum(base_mu[0:3])
    base_mu1 = base_mu[0]
    base_mu2 = base_mu[1]
    base_mu3 = base_mu[2]
    base_mu4 = base_mu[3]
    core_site_mu1 = 0.25
    core_site_mu2 = 0.25
    core_site_mu3 = 0.25
    core_site_mu4 = 0.25

    dist_mat = gen_distances_elfi(size_core, size_pan, core_mu, avg_gene_freq, prop_gene, gene_gl,
                                  base_mu1, base_mu2, base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3,
                                  core_site_mu4, pop_size, n_gen, max_hamming_core, max_jaccard_acc, True)

    dist_mat[:, 0] = dist_mat[:, 0] * max_hamming_core
    dist_mat[:, 1] = dist_mat[:, 1] * max_jaccard_acc

    return index, dist_mat

if __name__ == "__main__":
    # distfile = "distances/GPSv4_distances_sample1.txt"
    # threads = 1
    # paramsfile = "parameter_example.txt"
    # outpref = "test"

    options = get_options()
    threads = options.threads
    distfile = options.distfile
    paramsfile = options.params
    outpref = options.outpref

    # read in real files
    df = read_distfile(distfile)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(df["Core"].max())
    max_jaccard_acc = float(df["Accessory"].max())
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))

    params_list = []
    with open(paramsfile, "r") as f:
        next(f)
        for line in f:
            params_list.append(line.rstrip().split("\t"))

    print("Running Simulations")
    with Pool(processes=threads) as pool:
        for index, dist_mat in tqdm.tqdm(pool.imap(
                partial(run_sim, params_list=params_list, max_hamming_core=max_hamming_core, max_jaccard_acc=max_jaccard_acc),
                range(0, len(params_list))), total=len(params_list)):
            fig, ax = plt.subplots()

            x = df["Core"]
            y = df["Accessory"]

            ax.scatter(x, y, label="Real")

            x = np.array(dist_mat[:, 0])
            y = np.array(dist_mat[:, 1])

            ax.scatter(x, y, label="Sim " + str(index + 1))

            ax.set_xlabel("Core Hamming")
            ax.set_ylabel("Accessory Jaccard")

            fig.savefig(outpref + "_" + "sim_" + str(index + 1) + ".png")

            ax.clear()
            plt.clf
            plt.cla


