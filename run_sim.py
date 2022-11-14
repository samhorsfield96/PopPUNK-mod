import argparse
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
                    help='Size of core genome alignment (in bases). Default = 1140000 ')
    IO.add_argument('--core-var',
                    type=int,
                    default=106196,
                    help='Number of variant sites in core. Default = 106196 ')
    IO.add_argument('--base-freq',
                    default="0.25,0.25,0.25,0.25",
                    help='Base frequencies in starting core genome in order "A,C,G,T". Default = "0.25,0.25,0.25,0.25" ')
    IO.add_argument('--base-mu',
                    default="0.25,0.25,0.25,0.25",
                    help='Mutation rates from all other bases to each base, in order "A,C,G,T". Default = "0.25,0.25,0.25,0.25" ')
    IO.add_argument('--start-gene-freq',
                    default="0.5,0.5",
                    help='Gene frequencies in starting accessory genome in order "0,1". Default = "0.5,0.5" ')
    IO.add_argument('--avg-gene-freq',
                    type=float,
                    default=0.5,
                    help='Average gene frequency in accessory genome. Default = "0.5" ')
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
                    help='Dispersion for simulated core values as proportion of max mean value. '
                         'Default = 0.01')
    IO.add_argument('--sim-acc-dispersion',
                    default=0.01,
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
                    default=1,
                    help='Number of simulations to run. Default = 1')
    IO.add_argument('--adjust',
                    default=False,
                    action="store_true",
                    help='Adjust core and accessory distances for invariant sites. Default = False')
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
    size_core = options.core_size
    num_core = options.num_core
    num_pangenome = options.num_pan
    core_num_var = options.core_var
    adjusted = options.adjust
    avg_gene_freq = options.avg_gene_freq
    sim_core_dispersion = options.sim_core_dispersion
    sim_acc_dispersion = options.sim_acc_dispersion
    core_sites_man = options.core_sites_man
    acc_sites_man = options.acc_sites_man

    # determine per-site substitution rates
    core_sites_man = None
    acc_sites_man = None
    if options.core_sites_man is not None:
        core_sites_man = [float(i) for i in options.core_sites_man.split(",")]
    if options.acc_sites_man is not None:
        acc_sites_man = [float(i) for i in options.acc_sites_man.split(",")]

    # calculate number of core invariant sites
    core_num_invar = size_core - core_num_var

    # core mu is number of differences per base of alignment. Divide by two to account for using two diverging sequences
    core_mu = []
    str_core_mu = [float(i) / 2 for i in options.core_mu.split(",")]
    mu = str_core_mu[0]
    while mu <= str_core_mu[1]:
        core_mu.append(mu)
        mu += str_core_mu[2]

    # round to 6 dp
    core_mu = [round(i, 6) for i in core_mu]
    core_mu_ori = core_mu.copy()

    # to ensure dispersion is constant, set on max value of core
    sim_core_dispersion = max(core_mu) * sim_core_dispersion

    # adjust core mutation rate in variable sites to match overall core_mu
    #core_mu = [(i * (size_core / core_num_var)) for i in core_mu]

    print("Individual genome core per-base mutation rate values: ")
    print(core_mu)

    # acc mu is number of differences per gene in accessory genome. Divide by two to account for using two diverging sequences
    acc_mu = []
    str_acc_mu = [float(i) / 2 for i in options.acc_mu.split(",")]
    mu = str_acc_mu[0]
    while mu <= str_acc_mu[1]:
        acc_mu.append(mu)
        mu += str_acc_mu[2]

    #round to 6 dp
    acc_mu = [round(i, 6) for i in acc_mu]

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

    hamming_core_sims = [None] * options.num_sim
    hamming_acc_sims = [None] * options.num_sim
    jaccard_core_sims = [None] * options.num_sim
    jaccard_acc_sims = [None] * options.num_sim
    acc_vs_core_sims = [None] * options.num_sim
    pangenome_frac_sims = [None] * options.num_sim

    sim_list = []

    for i in range(options.num_sim):
        core_ref = np.random.choice(base_choices, size_core, p=base_freq)
        acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)

        ref_name = str(i) + "_0"

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

        hamming_core = [None] * len(core_mu)
        hamming_acc = [None] * len(acc_mu)
        jaccard_core = [None] * len(core_mu)
        jaccard_acc = [None] * len(acc_mu)
        acc_vs_core = [None] * len(acc_mu)
        pangenome_frac = [None] * len(acc_mu)

        print("Simulation " + str(i + 1))
        with Pool(processes=threads) as pool:
            for ind, hcore, hacc, jcore, jacc, avc, p_frac in tqdm.tqdm(pool.imap(
                    partial(gen_distances, core_var=core_var, acc_var=acc_ref, core_invar=core_invar,
                            num_core=num_core, core_mu=core_mu, acc_mu=acc_mu, adj=adjusted, avg_gene_freq=avg_gene_freq,
                            base_mu=base_mu, core_site_mu=core_site_mu, acc_site_mu=acc_site_mu,
                            sim_core_dispersion=sim_core_dispersion, sim_acc_dispersion=sim_acc_dispersion),
                    range(0, len(core_mu))), total=len(core_mu)):
                hamming_core[ind] = hcore
                hamming_acc[ind] = hacc
                jaccard_core[ind] = jcore
                jaccard_acc[ind] = jacc
                acc_vs_core[ind] = avc
                pangenome_frac[ind] = p_frac
                query_name = str(i) + "_" + str(ind + 1)
                sim_list.append((ref_name, query_name, hcore, jacc, avc))

        hamming_core_sims[i] = hamming_core
        hamming_acc_sims[i] = hamming_acc
        jaccard_core_sims[i] = jaccard_core
        jaccard_acc_sims[i] = jaccard_acc
        acc_vs_core_sims[i] = acc_vs_core
        pangenome_frac_sims[i] = pangenome_frac

    # print simulated run
    with open(options.outpref + "simulation.txt", "w") as f:
        for entry in sim_list:
            f.write(str(entry[0]) + "\t" + str(entry[1]) + "\t" + str(entry[2]) + "\t" + str(entry[3]) + "\t" + str(entry[4]) + "\n")

    print("Actual accessory : core rates")
    print(acc_vs_core_sims)

    print("Actual pangenome fractions")
    print(pangenome_frac_sims)

    print("Generating graphs...")

    mu_rates = (core_mu_ori, acc_mu_ori, core_mu_ori, acc_mu_ori)
    distances = (hamming_core_sims, hamming_acc_sims, jaccard_core_sims, jaccard_acc_sims)
    mu_names = ("core_mu", "acc_mu", "core_mu", "acc_mu")
    distance_names = ("hamming_core", "hamming_acc", "jaccard_core", "jaccard_acc")
    lengths = (size_core, num_pangenome, size_core, num_pangenome)

    # plot pangenome_fracs against accessory and core distances
    check_panfrac((distances[0], distances[3]), (pangenome_frac_sims, pangenome_frac_sims), options.outpref)

    # calculate proportion of invariable sites
    core_adj = size_core / core_num_var
    acc_adj = num_pangenome / size_acc

    generate_graph(mu_rates, distances, mu_names, distance_names, options.outpref, core_adj, 1, adjusted)

    print("Done.")