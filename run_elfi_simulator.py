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
    description = 'Run simulator of gene gain model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_elfi_simulator.py')

    IO = parser.add_argument_group('Input/Output options')
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

def run_sim(index, params_list, max_real_core, max_hamming_core, max_jaccard_acc):
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

    # set evenly spaced core hamming values across generations. Divide by two as pairwise divergence
    core_mu = (max_real_core / n_gen) / 2

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

    dist_mat, avg_core, avg_acc = gen_distances_elfi(size_core, size_pan, core_mu, avg_gene_freq, prop_gene, gene_gl,
                                                     base_mu1, base_mu2, base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3,
                                                     core_site_mu4, pop_size, n_gen, max_hamming_core, max_jaccard_acc, True)

    dist_mat[:, 0] = dist_mat[:, 0] * max_hamming_core
    dist_mat[:, 1] = dist_mat[:, 1] * max_jaccard_acc

    return index, dist_mat, avg_core, avg_acc

if __name__ == "__main__":
    # distfile = "distances/GPSv4_distances_sample1.txt"
    # threads = 1
    # paramsfile = "parameter_example_test.txt"
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
    #mode_hamming_core = float(df["Core"].mode())

    params_list = []
    with open(paramsfile, "r") as f:
        next(f)
        for line in f:
            params_list.append(line.rstrip().split("\t"))

    print("Running Simulations...")
    with Pool(processes=threads) as pool:
        for index, dist_mat, avg_core, avg_acc in tqdm.tqdm(pool.imap(
                partial(run_sim, params_list=params_list, max_real_core=max_real_core, max_hamming_core=max_hamming_core, max_jaccard_acc=max_jaccard_acc),
                range(0, len(params_list))), total=len(params_list)):
            # save to file
            np.savetxt(outpref + "_results_sim_" + str(index + 1) + ".csv", dist_mat, delimiter=",")
            avg_mat = np.zeros((avg_core.shape[0], 2))
            avg_mat[:, 0] = avg_core
            avg_mat[:, 1] = avg_acc
            np.savetxt(outpref + "_averages_sim_" + str(index + 1) + ".csv", avg_mat, delimiter=",")

            fig, ax = plt.subplots()

            x = df["Core"]
            y = df["Accessory"]

            ax.scatter(x, y, label="Real")

            x = np.array(dist_mat[:, 0])
            y = np.array(dist_mat[:, 1])

            ax.scatter(x, y, label="Sim " + str(index + 1))

            ax.set_xlabel("Core Hamming")
            ax.set_ylabel("Accessory Jaccard")

            fig.savefig(outpref + "_" + "sim_" + str(index + 1) + ".svg")

            ax.clear()
            plt.clf
            plt.cla

            # plot average distances over time
            y = avg_core

            ax.plot(y)

            ax.set_xlabel("Generation")
            ax.set_ylabel("Average Core Hamming")

            fig.savefig(outpref + "_" + "avg_core_sim_" + str(index + 1) + ".svg")

            ax.clear()
            plt.clf
            plt.cla

            # plot average distances over time
            y = avg_acc

            ax.plot(y)

            ax.set_xlabel("Generation")
            ax.set_ylabel("Average Acc Jaccard")

            fig.savefig(outpref + "_" + "avg_acc_sim_" + str(index + 1) + ".svg")

            ax.clear()
            plt.clf
            plt.cla


