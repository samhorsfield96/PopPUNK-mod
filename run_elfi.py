import argparse

import numpy as np
import scipy.stats
import logging
logging.basicConfig(level=logging.INFO)

import elfi
from simulate_divergence import *
import os
import sys

def get_options():
    description = 'Fit model to PopPUNK data using Approximate Baysesian computation'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_ELFI.py')

    IO = parser.add_argument_group('Input/Output options')
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
    IO.add_argument('--batch-size',
                    type=int,
                    default=10000,
                    help='Batch size for processing. Default = 10000 ')
    IO.add_argument('--samples',
                    type=int,
                    default=1000,
                    help='No. samples for posterior estimation. Default = 1000 ')
    IO.add_argument('--schedule',
                    type=str,
                    default="0.7,0.2,0.05",
                    help='SMC schedule, a list of thresholds to use for each population. Default = 0.7,0.2,0.05 ')
    IO.add_argument('--qnt',
                    type=float,
                    default=0.01,
                    help='Quantile of the samples with smallest discrepancies is accepted. Default = 0.01 ')
    IO.add_argument('--init-evidence',
                    type=int,
                    default=5000,
                    help='Number of initialization points sampled straight from the priors before starting to '
                         'optimize the acquisition of points. Default = 5000 ')
    IO.add_argument('--update-int',
                    type=int,
                    default=10,
                    help='Defines how often the GP hyperparameters are optimized. '
                         'Default = 10 ')
    IO.add_argument('--acq-noise-var',
                    type=float,
                    default=0.1,
                    help='Defines the diagonal covariance of noise added to the acquired points. '
                         'Default = 0.1 ')
    IO.add_argument('--n-evidence',
                    type=int,
                    default=5000,
                    help='Evidence points requested (including init-evidence). '
                         'Default = 5000 ')
    IO.add_argument('--distfile',
                    required=True,
                    help='popPUNK distance file to fit to. ')
    IO.add_argument('--seed',
                    type=int,
                    default=254,
                    help='Seed for random number generation. Default = 254. ')
    IO.add_argument('--summary',
                    choices=['quantile', 'mean'],
                    default="quantile",
                    help='Mode for summary statistics, either "mean" or "quantile". Default = "quantile". ')
    IO.add_argument('--mode',
                    choices=['rejection', 'SMC', 'BOLFI'],
                    default="rejection",
                    help='Mode for running model fit, either "rejection", "SMC or "BOLFI". Default = "rejection". ')
    IO.add_argument('--outpref',
                    default="PopPUNK-mod",
                    help='Output prefix. Default = "PopPUNK-mod"')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')
    IO.add_argument('--cluster',
                    action='store_true',
                    default=False,
                    help='Parallelise using ipyparallel if using cluster. Default = False')

    return parser.parse_args()

def mean(x):
    C = np.mean(x, axis=1)
    return C

def median(x):
    C = np.median(x, axis=1)
    return C

def max(x):
    C = np.max(x, axis=1)
    return C

def min(x):
    C = np.min(x, axis=1)
    return C

def quantile(x, q):
    C = np.quantile(x, q, axis=1)
    return C

def stddev(x):
    C = np.std(x, axis=1)
    return C

def gen_distances_elfi(size_core, size_pan, core_mu, avg_gene_freq, prop_gene, gene_gl,
                       base_mu1, base_mu2, base_mu3, base_mu4,
                       core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4,
                       pop_size, n_gen, batch_size=1, random_state=None):
    # determine vectors of core and accessory per-site mutation rate
    core_mu_arr = np.array([core_mu] * batch_size)
    acc_mu_arr = core_mu_arr * gene_gl

    # generate vectors for mutation rates
    base_mu = np.tile(np.array([base_mu1, base_mu2, base_mu3, base_mu4]), (batch_size, 1))
    core_site = np.tile(np.array([core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4]), (batch_size, 1))
    acc_site_1 = prop_gene
    acc_site_2 = 1 - acc_site_1
    acc_site = np.stack((acc_site_1, acc_site_2), axis=1)
    #print(acc_site)
    gene_mu = np.stack(([1 - avg_gene_freq] * batch_size, [avg_gene_freq] * batch_size), axis=1)
    #print(gene_mu)

    # simulate population forward using fisher-wright

    # TODO need a way to determine which sites will be fast or slow
    core_site_mu = calc_man_vec(size_core, size_core, core_site, batch_size)
    acc_site_mu = calc_man_vec(size_pan, size_pan, acc_site, batch_size)
    #print(core_site_mu)
    #print(acc_site_mu)

    # calculate accessory total gene gain/loss

    # generate starting genomes, rows are batches, columns are positions
    core_ref = np.zeros((batch_size, size_core))
    acc_ref = np.zeros((batch_size, size_pan))
    for i in range(batch_size):
        core_ref[i] = np.random.choice([1, 2, 3, 4], size_core, p=base_mu[i])
        acc_ref[i] = np.random.choice([0, 1], size_pan, p=gene_mu[i])

    pop_core = np.array([core_ref] * pop_size)
    pop_acc = np.array([acc_ref] * pop_size)

    # simulate population forward using fisher-wright
    for gen in range(1, n_gen):
        # sample from previous generation in each batch with replacement
        #print("before")
        #print(pop_core)
        if gen > 1:
            sample = np.random.choice(pop_core.shape[2], pop_size, replace=True)
            #print(sample)
            pop_core = pop_core[sample, :, :]
            pop_acc = pop_acc[sample, :, :]
        #print("after")
        #print(pop_core)

        # mutate genomes
        #print("before")
        #print(pop_core)
        pop_core = sim_divergence_vec(pop_core, core_mu_arr, True, base_mu, core_site_mu, pop_size)
        #print("after")
        #print(pop_core_new)
        pop_acc = sim_divergence_vec(pop_acc, acc_mu_arr, False, gene_mu, acc_site_mu, pop_size)


    for j in range(0, batch_size):
        pop_core_slice = pop_core[:, j, :]
        pop_acc_slice = pop_acc[:, j, :]
        #print(pop_core_slice)
        #print(pop_acc_slice)
        eucl = []

        # iterate over all genomes in population, calculating gamming distance
        for k in range(0, pop_size):
            for l in range(0, pop_size):
                if l < k:
                    hamming_core = distance.hamming(pop_core_slice[k], pop_core_slice[l])
                    jaccard_acc = distance.jaccard(pop_acc_slice[k], pop_acc_slice[l])
                    eucl.append(math.sqrt((hamming_core ** 2) + (jaccard_acc ** 2)))

        #print(eucl)
        if j == 0:
            eucl_mat = np.zeros((batch_size, len(eucl)))
        eucl_mat[j] = np.array(eucl)

    return eucl_mat

if __name__ == "__main__":
    #testing
    # size_core = 2
    # size_pan = 3
    # avg_gene_freq = 0.5
    # batch_size = 10
    # N_samples = 10
    # qnt = 0.01
    # seed = 254
    # summary = "quantile"
    # data_dir = "distances"
    # data_pref = "Pneumo_sim_simulation"
    # num_steps = 10
    # threads = 4
    # mode = "rejection"
    # outpref = "test"
    # initial_evidence = 20
    # update_interval = 10
    # acq_noise_var = 0.1
    # n_evidence = 200
    # info_freq = 1000
    # base_mu = [0.25, 0.25, 0.25, 0.25]
    # cluster = False
    # complexity = "simple"
    # schedule = "0.7,0.2,0.05"
    # pop_size = 5
    # n_gen = 100

    options = get_options()
    threads = options.threads
    distfile = options.distfile
    size_core = options.core_size
    size_pan = options.pan_size
    batch_size = options.batch_size
    qnt = options.qnt
    N_samples = options.samples
    seed = options.seed
    outpref = options.outpref
    summary = options.summary
    mode = options.mode
    initial_evidence = options.init_evidence
    update_interval = options.update_int
    acq_noise_var = options.acq_noise_var
    n_evidence = options.n_evidence
    avg_gene_freq = options.avg_gene_freq
    base_mu = [float(i) for i in options.base_mu.split(",")]
    cluster = options.cluster
    schedule = options.schedule
    n_gen = options.ngen
    pop_size = options.pop_size

    # parse schedule
    schedule = [float(x) for x in schedule.split(",")]

    # read in real files
    df = read_file(distfile)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(df["Core"].max())
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))

    # set constants
    # set evenly spaced core hamming values across generations
    core_mu = max_real_core / n_gen

    # set minimum prop_core_var and prop_acc_var based on number of sequence bins (hard coded at 4 at the moment)
    min_prop_core_var = 4 / size_core
    min_prop_acc_var = 2 / size_pan

    # set priors
    # priors for gene gain and loss rates per site
    max_value = 10**6
    gene_gl = elfi.Prior('uniform', 0, max_value)

    # prior for size two gene compartments
    prop_gene = elfi.Prior('uniform', 0.5, 1 - 0.5)

    # round to 6 dp
    base_mu = [round(i, 6) for i in base_mu]

    # ensure probabilities sum to 1
    if sum(base_mu) != 1:
        base_mu[-1] = 1 - sum(base_mu[0:3])
    base_mu1 = base_mu[0]
    base_mu2 = base_mu[1]
    base_mu3 = base_mu[2]
    base_mu4 = base_mu[3]
    # base_mu1 = elfi.Prior('uniform', 0, 1)
    # base_mu2 = elfi.Prior('uniform', 0, 1)
    # base_mu3 = elfi.Prior('uniform', 0, 1)
    # base_mu4 = elfi.Prior('uniform', 0, 1)

    core_site_mu1 = 0.25
    core_site_mu2 = 0.25
    core_site_mu3 = 0.25
    core_site_mu4 = 0.25
    # core_site_mu1 = elfi.Prior('uniform', 0, 1)
    # core_site_mu2 = elfi.Prior('uniform', 0, 1)
    # core_site_mu3 = elfi.Prior('uniform', 0, 1)
    # core_site_mu4 = elfi.Prior('uniform', 0, 1)
    #core_site_mu5 = elfi.Prior('uniform', 0, 1)

    #get observed data
    obs_core = df['Core'].to_numpy()
    obs_acc = df['Accessory'].to_numpy()

    # calculate euclidean distance to origin
    obs = np.sqrt((obs_core ** 2) + (obs_acc ** 2)).reshape(1, -1)

    Y = elfi.Simulator(gen_distances_elfi, size_core, size_pan, core_mu, avg_gene_freq, prop_gene, gene_gl,
                       base_mu1, base_mu2, base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3,
                       core_site_mu4, pop_size, n_gen, observed=obs)

    S_min = elfi.Summary(min, Y)
    S_max = elfi.Summary(max, Y)
    if summary == "quantile":
        # generate summary statitics as quantiles of data
        S_q1 = elfi.Summary(quantile, Y, 0.1)
        S_q2 = elfi.Summary(quantile, Y, 0.2)
        S_q3 = elfi.Summary(quantile, Y, 0.3)
        S_q4 = elfi.Summary(quantile, Y, 0.4)
        S_q5 = elfi.Summary(quantile, Y, 0.5)
        S_q6 = elfi.Summary(quantile, Y, 0.6)
        S_q7 = elfi.Summary(quantile, Y, 0.7)
        S_q8 = elfi.Summary(quantile, Y, 0.8)
        S_q9 = elfi.Summary(quantile, Y, 0.9)

        d = elfi.Distance('euclidean', S_min, S_q1, S_q2, S_q3, S_q4, S_q5, S_q6, S_q7, S_q8, S_q9, S_max)
    else:
        S_mean = elfi.Summary(mean, Y)
        S_stddev = elfi.Summary(stddev, Y)
        d = elfi.Distance('euclidean', S_min, S_mean, S_stddev, S_max)

    #set multiprocessing client
    if cluster == True:
        # must start ipyarallel cluster e.g. !ipcluster start -n threads --daemon
        elfi.set_client('ipyparallel')
    else:
        if threads > 1:
            elfi.set_client('multiprocessing')
            elfi.set_client(elfi.clients.multiprocessing.Client(num_processes=threads))
        else:
            elfi.set_client('native')

    os.environ['NUMEXPR_MAX_THREADS'] = str(threads)

    if mode == "rejection":
        rej = elfi.Rejection(d, batch_size=batch_size, seed=seed)

        result = rej.sample(N_samples, quantile=qnt)
    elif mode == "SMC":
        smc = elfi.SMC(d, batch_size=batch_size, seed=seed)
        result = smc.sample(N_samples, schedule)
    else:
        log_d = elfi.Operation(np.log, d)
        bounds = {
            'gene_gl' : (0, max_value),
            'prop_gene' : (0, 1),
        }
        bolfi = elfi.BOLFI(log_d, batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                           acq_noise_var=acq_noise_var, seed=seed, bounds=bounds)
        post = bolfi.fit(n_evidence=n_evidence)
        result = bolfi.sample(N_samples)

    with open(outpref + "_ELFI_summary.txt", "w") as f:
        print(result, file=f)

    # plot graphs
    # plot marginals
    result.plot_marginals()
    plt.savefig(outpref + '_marginals.png')

    plt.clf
    plt.cla

    # plot paired marginals
    result.plot_pairs()
    plt.savefig(outpref + '_pairs.png')

    plt.clf
    plt.cla

    # plot BOLFI specific graphs
    if mode == "BOLFI":
        #plot MCMC traces
        result.plot_traces();
        plt.savefig(outpref + '_BOLFI_traces.png')

        plt.clf
        plt.cla

        # plot discrepancies
        bolfi.plot_discrepancy()

        plt.savefig(outpref + '_BOLFI_discepancy.png')

        plt.clf
        plt.cla
