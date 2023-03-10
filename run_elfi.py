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
    IO.add_argument('--load',
                    default=None,
                    help='Directory of previous ELFI model. Required if running "sample" mode ')
    IO.add_argument('--seed',
                    type=int,
                    default=254,
                    help='Seed for random number generation. Default = 254. ')
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


def gen_distances_elfi(size_core, size_pan, core_mu, avg_gene_freq, ratio_gene_gl, gene_gl_speed, prop_gene,
                       base_mu1, base_mu2, base_mu3, base_mu4,
                       core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4,
                       pop_size, n_gen, max_real_core, simulate, batch_size=1, random_state=None):
    # determine vectors of core and accessory per-site mutation rate.
    core_mu_arr = np.array([core_mu] * batch_size)
    acc_mu_arr = core_mu_arr * gene_gl_speed

    # core mu array increased by factor max_real_core as only looking at subset
    #core_mu_arr = core_mu_arr / max_real_core

    # calculate actual number of sites mutating
    #size_core_mut = round(max_real_core * size_core)
    size_core_mut = size_core

    # generate vectors for mutation rates
    base_mu = np.tile(np.array([base_mu1, base_mu2, base_mu3, base_mu4]), (batch_size, 1))
    core_site = np.tile(np.array([core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4]), (batch_size, 1))
    acc_site_1 = ratio_gene_gl
    acc_site_2 = 1 - acc_site_1
    if simulate:
        acc_site = np.stack(([acc_site_1], [acc_site_2]), axis=1)
        proportion_gene = np.stack(([prop_gene], [1 - prop_gene]), axis=1)
    else:
        acc_site = np.stack((acc_site_1, acc_site_2), axis=1)
        proportion_gene = np.stack((prop_gene, 1 - prop_gene), axis=1)
    gene_mu = np.stack(([1 - avg_gene_freq] * batch_size, [avg_gene_freq] * batch_size), axis=1)

    # calculate per-site mutation rate
    core_site_mu = calc_man_vec(size_core_mut, size_core_mut, core_site, batch_size)
    acc_site_mu = calc_man_vec(size_pan, size_pan, acc_site, batch_size, proportion_gene)

    # generate starting genomes, rows are batches, columns are positions
    core_ref = np.zeros((batch_size, size_core_mut))
    acc_ref = np.zeros((batch_size, size_pan))
    for i in range(batch_size):
        core_ref[i] = np.random.choice([1, 2, 3, 4], size_core_mut, p=base_mu[i])
        acc_ref[i] = np.random.choice([0, 1], size_pan, p=gene_mu[i])

    pop_core = np.repeat([core_ref.copy()], pop_size, axis=0)
    pop_acc = np.repeat([acc_ref.copy()], pop_size, axis=0)

    # generate core tuple
    choices_1 = np.array([2, 3, 4])
    choices_2 = np.array([1, 3, 4])
    choices_3 = np.array([1, 2, 4])
    choices_4 = np.array([1, 2, 3])
    prob_1 = base_mu.copy()
    prob_2 = base_mu.copy()
    prob_3 = base_mu.copy()
    prob_4 = base_mu.copy()
    prob_1 = np.delete(prob_1, 0, 1)
    prob_1 = prob_1 / np.sum(prob_1)
    prob_2 = np.delete(prob_2, 1, 1)
    prob_2 = prob_2 / np.sum(prob_2)
    prob_3 = np.delete(prob_3, 2, 1)
    prob_3 = prob_3 / np.sum(prob_3)
    prob_4 = np.delete(prob_4, 3, 1)
    prob_4 = prob_4 / np.sum(prob_4)

    core_tuple = (choices_1, choices_2, choices_3, choices_4, prob_1, prob_2, prob_3, prob_4)

    # run numba-backed WF model
    pop_core, pop_acc, avg_core, avg_acc = run_WF_model(pop_core, pop_acc, n_gen, pop_size, core_mu_arr, acc_mu_arr,
                                                        core_site_mu, acc_site_mu, max_real_core, simulate, core_tuple)

    # run numba-backed distance calculator
    core_mat, acc_mat = calc_dists(pop_core, pop_acc, batch_size, max_real_core, simulate)

    if simulate:
        dist_mat = np.zeros((core_mat.shape[0], 2))
        dist_mat[:, 0] = core_mat
        dist_mat[:, 1] = acc_mat
        return dist_mat, avg_core, avg_acc
    else:
        dist_mat = np.zeros((batch_size, (acc_mat.shape[1])))
        for j in range(0, batch_size):
            #dist_mat[j] = np.concatenate([core_mat[j], acc_mat[j]])
            dist_mat[j] = acc_mat
        return dist_mat

if __name__ == "__main__":
    #testing
    # size_core = 4
    # size_pan = 2
    # avg_gene_freq = 0.5
    # batch_size = 10
    # N_samples = 10
    # qnt = 0.01
    # seed = 254
    # distfile = "distances/GPSv4_distances_sample1.txt"
    # num_steps = 10
    # threads = 4
    # mode = "BOLFI"
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
    # load = "test_pools/outputpool_254"
    # run_mode = "sim"

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
    load = options.load
    run_mode = options.run_mode

    # parse schedule
    schedule = [float(x) for x in schedule.split(",")]

    # read in real files
    df = read_distfile(distfile)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(df["Core"].max())
    max_jaccard_acc = float(df["Accessory"].max())
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))

    # set constants
    # set evenly spaced core hamming values across generations
    core_mu = (max_real_core / (n_gen - 1)) / 2

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

    #get observed data, normalise
    #obs_core = get_quantile(df['Core'].to_numpy())# / max_hamming_core)
    #obs_acc = get_quantile(df['Accessory'].to_numpy())# / max_jaccard_acc)
    obs_acc = np.histogram(df['Accessory'].to_numpy(), bins=50, range=(0, 1), density=True)[0]

    # calculate euclidean distance to origin
    #obs = np.concatenate([obs_core, obs_acc])
    obs = obs_acc

    # set priors
    # priors for gene gain and loss rates per site
    max_value = 10 ** 6
    gene_gl_speed = elfi.Prior('uniform', 0, max_value)

    # prior for difference in probability of sample fast vs. slow two gene compartments
    ratio_gene_gl = elfi.Prior('uniform', 0.5, 1 - 0.5)

    # prior for difference size of compartments
    prop_gene = elfi.Prior('uniform', 0, 1)

    Y = elfi.Simulator(gen_distances_elfi, size_core, size_pan, core_mu, avg_gene_freq, ratio_gene_gl, gene_gl_speed, prop_gene,
                        base_mu1, base_mu2, base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3,
                        core_site_mu4, pop_size, n_gen, max_real_core, False, observed=obs)

    d = elfi.Distance('euclidean', Y)

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

    if run_mode == "sim":
        print("Simulating data...")
        if mode == "rejection":
            mod = elfi.Rejection(d, batch_size=batch_size, seed=seed)
            result = mod.sample(N_samples, quantile=qnt)

        elif mode == "SMC":
            mod = elfi.SMC(d, batch_size=batch_size, seed=seed)
            result = mod.sample(N_samples, schedule)
        else:
            log_d = elfi.Operation(np.log, d)
            bounds = {
                'gene_gl_speed' : (0, max_value),
                'ratio_gene_gl' : (0.5, 1),
                'prop_gene' : (0, 1),
            }
            mod = elfi.BOLFI(log_d, batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                               acq_noise_var=acq_noise_var, seed=seed, bounds=bounds)

            post = mod.fit(n_evidence=n_evidence)
            result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence)

            # not implemented for more than 2 dimensions
            # post.plot(logpdf=True)
            # plt.savefig("posterior.svg")
            # plt.close()

            mod.plot_discrepancy()
            plt.savefig(outpref + "_BOLFI_discrepancy.svg")
            plt.close()

            # plot MCMC traces
            result.plot_traces();
            plt.savefig(outpref + '_BOLFI_traces.svg')
            plt.close()

            # # plot results
            # mod.plot_state()
            # plt.savefig(outpref + "_" + mode + "_state.svg")
            # plt.close()

        with open(outpref + "_ELFI_summary.txt", "w") as f:
            print(result, file=f)

        # save model
        save_path = outpref + '_pools'
        os.makedirs(save_path, exist_ok=True)
        arraypool = elfi.OutputPool(['ratio_gene_gl', 'gene_gl_speed', 'prop_gene', 'Y', 'd'], prefix=save_path)
        arraypool.set_context(mod)
        arraypool.save()

    else:
        print("Loading models in {}".format(load))
        if load == None:
            print('Previously saved ELFI pool required for "sample" mode. Please specify "--load."')
            sys.exit(1)

        # parse filename
        load_pref = load.rsplit('/', 1)[0]
        load_name = load.rsplit('/', 1)[-1]

        arraypool = elfi.OutputPool.open(name=load_name, prefix=load_pref)

        if mode == "rejection":
            mod = elfi.Rejection(d, batch_size=batch_size, seed=seed, pool=arraypool)
            result = mod.sample(N_samples, quantile=qnt)
        elif mode == "SMC":
            mod = elfi.SMC(d, batch_size=batch_size, seed=seed, pool=arraypool)
            result = mod.sample(N_samples, schedule)
        else:
            log_d = elfi.Operation(np.log, d)
            bounds = {
                'gene_gl_speed' : (0, max_value),
                'ratio_gene_gl' : (0.5, 1),
                'prop_gene': (0, 1),
            }
            mod = elfi.BOLFI(log_d, batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                             acq_noise_var=acq_noise_var, seed=seed, bounds=bounds, pool=arraypool)

            result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence)

            mod.plot_discrepancy()
            plt.savefig(outpref + "_BOLFI_discrepancy.svg")
            plt.close()

            # plot results
            mod.plot_state()
            plt.savefig(outpref + "_" + mode + "_state.svg")
            plt.close()

            #plot MCMC traces
            result.plot_traces();
            plt.savefig(outpref + '_BOLFI_traces.svg')
            plt.close()

    with open(outpref + "_ELFI_summary.txt", "w") as f:
        print(result, file=f)

    # plot graphs
    # plot marginals
    result.plot_marginals()
    plt.savefig(outpref + "_" + mode + '_marginals.svg')

    plt.clf
    plt.cla

    # plot paired marginals
    result.plot_pairs()
    plt.savefig(outpref + "_" + mode + '_pairs.svg')
    plt.close()

    sys.exit(0)