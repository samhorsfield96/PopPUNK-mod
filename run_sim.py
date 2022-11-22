import argparse
import sys
from multiprocessing import Pool
from functools import partial
from simulate_divergence import *
import tqdm

def get_options():
    description = 'Calculate relationship between Hamming/Jaccard distances and core/accessory divergence'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_sim.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--core-size',
                    type=int,
                    default=1140000,
                    help='Size of core genome alignment (in bases). '
                         'Default = 1140000 ')
    IO.add_argument('--core-var',
                    type=int,
                    default=106196,
                    help='Number of variant sites in core. '
                         'Default = 106196 ')
    IO.add_argument('--base-freq',
                    default="0.3,0.2,0.2,0.3",
                    help='Base frequencies in starting core genome in order "A,C,G,T". '
                         'Default = "0.3,0.2,0.2,0.3" ')
    IO.add_argument('--base-mu',
                    default="0.3,0.2,0.2,0.3",
                    help='Mutation rates from all other bases to each base, in order "A,C,G,T". '
                         'Default = "0.3,0.2,0.2,0.3" ')
    IO.add_argument('--start-gene-freq',
                    default="0.608,0.392",
                    help='Gene frequencies in starting accessory genome in order "0,1". Default = "0.608,0.392" ')
    IO.add_argument('--avg-gene-freq',
                    type=float,
                    default=0.392,
                    help='Average gene frequency in accessory genome. Default = "0.392" ')
    IO.add_argument('--num-core',
                    type=int,
                    default=1194,
                    help='Number of core genes. Default = 1194 ')
    IO.add_argument('--num-pan',
                    type=int,
                    default=5442,
                    help='Number of genes in pangenome. Default = 5442')
    IO.add_argument('--core-mu',
                    default="0,0.02,0.002",
                    help='Range of core genome mutation rate values (mutations per site per genome) in form start,stop,step. '
                         'Default = "0,0.02,0.002"')
    IO.add_argument('--acc-mu',
                    default="0,0.4,0.04",
                    help='Range of accessory gene gain/loss rates (change per gene per genome) in form start,stop,step. '
                         'Default = "0,0.4,0.04"')
    IO.add_argument('--acc_func',
                    default=None,
                    help='Function for acc. variation around core. Takes form quadratic,slope,intercept. '
                         'Default = None')
    IO.add_argument('--core-sites',
                    type=int,
                    default=3,
                    help='Number of different core site mutation rates. '
                         'Default = 3')
    IO.add_argument('--acc-sites',
                    type=int,
                    default=3,
                    help='Number of different accessory site mutation rates. '
                         'Default = 3')
    IO.add_argument('--core-gamma-shape',
                    type=float,
                    default=20.0,
                    help='Shape parameter for core per-site substitution rates. '
                         'Default = 20.0')
    IO.add_argument('--core-gamma-scale',
                    type=float,
                    default=1.0,
                    help='Scale parameter for core per-site substitution rates. '
                         'Default = 1.0')
    IO.add_argument('--acc-gamma-shape',
                    type=float,
                    default=20.0,
                    help='Shape parameter for accessory per-site substitution rates. '
                         'Default = 20.0')
    IO.add_argument('--acc-gamma-scale',
                    type=float,
                    default=1.0,
                    help='Scale parameter for accessory per-site substitution rates. '
                         'Default = 1.0')
    IO.add_argument('--sim-core-dispersion',
                    default=0.01,
                    type=float,
                    help='Dispersion for simulated core values as proportion of max mean value. '
                         'Default = 0.01')
    IO.add_argument('--sim-acc-dispersion',
                    default=0.01,
                    type=float,
                    help='Dispersion for simulated accessory values as proportion of max mean value. '
                         'Default = 0.01')
    IO.add_argument('--core-sites-man',
                    default=None,
                    help='Manual core per-site mutation rates. Must sum to 1. '
                         'Default = None')
    IO.add_argument('--acc-sites-man',
                    default=None,
                    help='Manual accessory per-site mutation rates. Must sum to 1. '
                         'Default = None')
    IO.add_argument('--num-sim',
                    type=int,
                    default=5,
                    help='Number of simulations to run. Default = 5')
    IO.add_argument('--no-adjust',
                    default=True,
                    action="store_false",
                    help="Don't adjust core and accessory distances for invariant sites. ")
    IO.add_argument('--outpref',
                    default="popPUNK-mod",
                    help='Output prefix. Default = "popPUNK-mod"')
    IO.add_argument('--graph',
                    default=False,
                    action="store_true",
                    help="Generate distance comparison graphs. ")
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()
    threads = options.threads
    size_core = options.core_size
    num_core = options.num_core
    num_pangenome = options.num_pan
    core_num_var = options.core_var
    adjusted = options.no_adjust
    avg_gene_freq = options.avg_gene_freq
    sim_core_dispersion = options.sim_core_dispersion
    sim_acc_dispersion = options.sim_acc_dispersion
    core_sites_man = options.core_sites_man
    acc_sites_man = options.acc_sites_man
    acc_func = options.acc_func
    gen_graph = options.graph

    # determine per-site substitution rates
    core_sites_man = None
    acc_sites_man = None
    if options.core_sites_man is not None:
        core_sites_man = [float(i) for i in options.core_sites_man.split(",")]
    if options.acc_sites_man is not None:
        acc_sites_man = [float(i) for i in options.acc_sites_man.split(",")]

    # calculate number of core invariant sites
    core_num_invar = size_core - core_num_var

    # core mu is number of differences per base of alignment.
    core_mu = []
    str_core_mu = [float(i) for i in options.core_mu.split(",")]
    num_items = int(round((str_core_mu[1] / str_core_mu[2]), 0))
    core_mu.append(str_core_mu[0])
    for i in range(1, num_items):
        core_mu.append(i * str_core_mu[2])

    # round to 6 dp
    core_mu = [round(i, 6) for i in core_mu]
    core_mu_ori = core_mu.copy()

    # to ensure dispersion is constant, set on max value of core
    sim_core_dispersion = max(core_mu) * sim_core_dispersion

    # adjust core mutation rate in variable sites to match overall core_mu
    #core_mu = [(i * (size_core / core_num_var)) for i in core_mu]

    print("Individual genome core per-base mutation rate values: ")
    print(core_mu)

    # acc mu is number of differences per gene in accessory genome.
    acc_mu = []
    if acc_func is None:
        str_acc_mu = [float(i) for i in options.acc_mu.split(",")]
        num_items = int(round((str_acc_mu[1] / str_acc_mu[2]), 0))
        acc_mu.append(str_acc_mu[0])
        for i in range(1, num_items):
            acc_mu.append(i * str_acc_mu[2])
    else:
        str_acc_mu = [float(i) for i in options.acc_func.split(",")]
        acc_mu = [(str_acc_mu[0] * (mu ** 2)) + (mu * str_acc_mu[1]) + str_acc_mu[2] for mu in core_mu]

    #round to 6 dp
    acc_mu = [round(i, 6) for i in acc_mu]

    if any([x < 0  for x in acc_mu]):
        print("Parameters generated negative values")
        sys.exit(0)

    # to ensure dispersion is constant, set on max value of accessory
    sim_acc_dispersion = max(acc_mu) * sim_acc_dispersion

    # get number of accessory genes
    size_acc = num_pangenome - num_core
    acc_mu_ori = acc_mu.copy()

    # adjust acc_mu by factor of total pangenome / accessory pangenome size.
    #acc_mu = [(i * (num_pangenome / size_acc)) for i in acc_mu]

    print("Individual genome accessory per-gene mutation rate values: ")
    print(acc_mu)

    # parse base frequencies and generate core_ref
    base_freq = [float(i) for i in options.base_freq.split(",")]
    base_choices = [1, 2, 3, 4]

    # round to 6 dp
    base_freq = [round(i, 6) for i in base_freq]

    # ensure probabilities sum to 1
    if sum(base_freq) != 1:
        base_freq[-1] = 1 - sum(base_freq[0:3])

    print("Core base frequencies [A, C, G, T]: ")
    print(base_freq)

    # parse individual base mutation rates
    base_mu = [float(i) for i in options.base_mu.split(",")]

    # round to 6 dp
    base_mu = [round(i, 6) for i in base_mu]

    # ensure probabilities sum to 1
    if sum(base_mu) != 1:
        base_mu[-1] = 1 - sum(base_mu[0:3])

    # parse gene frequencies and generate acc_ref
    gene_freq = [float(i) for i in options.start_gene_freq.split(",")]
    gene_choices = [0, 1]

    # round to 6 dp
    gene_freq = [round(i, 6) for i in gene_freq]

    # ensure probabilities sum to 1
    if sum(gene_freq) != 1:
        base_freq[-1] = 1 - base_freq[0]

    print("Accessory gene frequencies [0, 1]: ")
    print(gene_freq)

    print("Iterating over distance simulations...")

    sim_list = []

    # start from same root genome
    core_ref = np.random.choice(base_choices, size_core, p=base_freq)
    acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)

    # pull out variable sites in core_ref
    sites = np.random.choice(range(core_ref.size), core_num_var, replace=False)
    present = np.full(core_ref.size, False)
    present[sites] = True
    core_var = core_ref[present]
    core_invar = core_ref[np.invert(present)]
    # determine per-site mutation rate
    if core_sites_man is None:
        if i == 0:
            print("Core genome per-site mutation rate compartments: ")
        core_site_mu = calc_gamma(core_var.size, options.core_sites, options.core_gamma_shape,
                                  options.core_gamma_scale, i)
    else:
        core_site_mu = calc_man(core_var.size, core_sites_man)

    if acc_sites_man is None:
        if i == 0:
            print("Accessory genome per-site mutation rate compartments: ")
        acc_site_mu = calc_gamma(size_acc, options.acc_sites, options.acc_gamma_shape,
                                 options.acc_gamma_scale, i)
    else:
        acc_site_mu = calc_man(size_acc, acc_sites_man)

    core_list = [None] * len(core_mu) * options.num_sim
    acc_list = [None] * len(acc_mu) * options.num_sim

    for i in range(options.num_sim):
        print("Simulation " + str(i + 1))
        with Pool(processes=threads) as pool:
            for ind, core_seq, acc_seq in tqdm.tqdm(pool.imap(
                    partial(gen_distances, core_var=core_var, acc_var=acc_ref, core_invar=core_invar,
                            num_core=num_core, core_mu=core_mu, acc_mu=acc_mu, adj=adjusted, avg_gene_freq=avg_gene_freq,
                            base_mu=base_mu, core_site_mu=core_site_mu, acc_site_mu=acc_site_mu,
                            sim_core_dispersion=sim_core_dispersion, sim_acc_dispersion=sim_acc_dispersion),
                    range(0, len(core_mu))), total=len(core_mu)):
                core_list[(len(core_mu) * i) + ind] = core_seq
                acc_list[(len(core_mu) * i) + ind] = acc_seq

    # go through all combinations and determine distances
    core_hamming = []
    core_jaccard = []
    acc_jaccard = []
    acc_hamming = []
    for i in range(0, len(core_list)):
        for j in range(0, len(core_list)):
            #ensure only part of array is searched
            if j <= i:
                continue
            core1 = core_list[i]
            core2 = core_list[j]
            acc1 = acc_list[i]
            acc2 = acc_list[j]

            hcore = distance.hamming(core1, core2)
            jacc = distance.jaccard(acc1, acc2)
            jcore = distance.jaccard(core1, core2)
            hacc = distance.hamming(acc1, acc2)

            core_hamming.append(hcore)
            acc_jaccard.append(jacc)
            acc_hamming.append(hacc)
            core_jaccard.append(jcore)

            sim_list.append((str(i), str(j), hcore, jacc))


    # print simulated run
    with open(options.outpref + "_simulation.txt", "w") as f:
        for entry in sim_list:
            f.write(str(entry[0]) + "\t" + str(entry[1]) + "\t" + str(entry[2]) + "\t" + str(entry[3]) + "\n")

    # print("Actual accessory : core rates")
    # print(acc_vs_core_sims)

    # print("Actual pangenome fractions")
    # print(pangenome_frac_sims)

    print("Generating graphs...")

    mu_rates = (core_mu_ori, acc_mu_ori, core_mu_ori, acc_mu_ori)
    distances = (core_hamming, acc_hamming, core_jaccard, acc_jaccard)
    mu_names = ("core_mu", "acc_mu", "core_mu", "acc_mu")
    distance_names = ("hamming_core", "hamming_acc", "jaccard_core", "jaccard_acc")
    lengths = (size_core, num_pangenome, size_core, num_pangenome)

    # plot pangenome_fracs against accessory and core distances
    #check_panfrac((distances[0], distances[3]), (pangenome_frac_sims, pangenome_frac_sims), options.outpref)

    # calculate proportion of invariable sites
    core_adj = size_core / core_num_var
    acc_adj = num_pangenome / size_acc

    generate_graph(mu_rates, distances, mu_names, distance_names, options.outpref, core_adj, 1, adjusted, gen_graph)

    print("Done.")
    sys.exit(0)