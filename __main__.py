import argparse
from multiprocessing import Pool
from functools import partial
from simulate_divergence import *
import tqdm

def get_options():
    description = 'Calculate relationship between Hamming/Jaccard distances and core/accessory divergence'
    parser = argparse.ArgumentParser(description=description,
                                     prog='distance_sim')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--aln-core',
                    type=int,
                    default=1140000,
                    help='Size of core genome alignment (in bases). Default = 1140000 ')
    IO.add_argument('--num-core',
                    type=int,
                    default=1194,
                    help='Number of core genes. Default = 1194 ')
    IO.add_argument('--num-pan',
                    type=int,
                    default=5442,
                    help='Number of genes in pangenome. Default = 5442')
    IO.add_argument('--core-mu',
                    default="0,2,0.2",
                    help='Range of core genome mutation rate values (mutations per site per genome) in form start,stop,step. '
                         'Default = "0,2,0.2"')
    IO.add_argument('--acc-mu',
                    default="0,2,0.2",
                    help='Range of accessory gene gain/loss rates (change per gene per genome) in form start,stop,step. '
                         'Default = "0,2,0.2"')
    IO.add_argument('--outpref',
                    default="./",
                    help='Output prefix. Default = "./"')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()
    threads = options.threads
    size_core = options.aln_core
    num_core = options.num_core
    num_pangenome = options.num_pan

    # core mu is number of differences per base of alignment
    core_mu = []
    str_core_mu = [float(i) for i in options.core_mu.split(",")]
    mu = str_core_mu[0]
    while mu <= str_core_mu[1]:
        core_mu.append(mu)
        mu += str_core_mu[2]

    # round to 6 dp
    core_mu = [round(i, 6) for i in core_mu]

    print("Core per-base mutation rate values: ")
    print(core_mu)

    # acc mu is number of differences per gene in accessory genome
    acc_mu = []
    str_acc_mu = [float(i) for i in options.acc_mu.split(",")]
    mu = str_acc_mu[0]
    while mu <= str_acc_mu[1]:
        acc_mu.append(mu)
        mu += str_acc_mu[2]

    #round to 6 dp
    acc_mu = [round(i, 6) for i in acc_mu]

    print("Accessory per-gene mutation rate values: ")
    print(acc_mu)

    hamming_core = [None] * len(core_mu)
    hamming_acc = [None] * len(acc_mu)
    jaccard_core = [None] * len(core_mu)
    jaccard_acc = [None] * len(acc_mu)

    print("Iterating over distance simulations...")

    with Pool(processes=threads) as pool:
        for ind, hcore, hacc, jcore, jacc in tqdm.tqdm(pool.imap(
            partial(gen_distances, size_core=size_core, num_core=num_core,
                    num_pangenome=num_pangenome, core_mu=core_mu, acc_mu=acc_mu),
                range(0, len(core_mu))), total = len(core_mu)):
            hamming_core[ind] = hcore
            hamming_acc[ind] = hacc
            jaccard_core[ind] = jcore
            jaccard_acc[ind] = jacc

    mu_rates = (core_mu, acc_mu, core_mu, acc_mu)
    distances = (hamming_core, hamming_acc, jaccard_core, jaccard_acc)
    mu_names = ("core_mu", "acc_mu", "core_mu", "acc_mu")
    distance_names = ("hamming_core", "hamming_acc", "jaccard_core", "jaccard_acc")

    print("Generating graphs...")

    generate_graph(mu_rates, distances, mu_names, distance_names, options.outpref)

    print("Done.")