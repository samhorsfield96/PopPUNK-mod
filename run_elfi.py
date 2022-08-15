import argparse
import scipy.stats

import elfi
from simulate_divergence import *
from fit_distances import read_files

def get_options():
    description = 'Calculate relationship between Hamming/Jaccard distances and core/accessory divergence'
    parser = argparse.ArgumentParser(description=description,
                                     prog='distance_sim')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--core-size',
                    type=int,
                    default=10000,
                    help='Number of positions in core genome. Default = 10000 ')
    IO.add_argument('--pan-size',
                    type=int,
                    default=1000,
                    help='Number of positions in pangenome. Default = 10000 ')
    IO.add_argument('--max-acc-vs-core',
                    type=int,
                    default=1000,
                    help='Maximum ratio between accessory and core genome evolution. Default = 1000 ')
    IO.add_argument('--num-steps',
                    type=int,
                    default=50,
                    help='Number of steps to take in increasing divergence. Default = 50 ')
    IO.add_argument('--batch-size',
                    type=int,
                    default=10000,
                    help='Batch size for processing. Default = 10000 ')
    IO.add_argument('--samples',
                    type=int,
                    default=1000,
                    help='No. samples for posterior estimation. Default = 1000 ')
    IO.add_argument('--qnt',
                    type=float,
                    default=0.01,
                    help='Quantile of the samples with smallest discrepancies is accepted. Default = 0.01 ')
    IO.add_argument('--init-evidence',
                    type=int,
                    default=262145,
                    help='Number of initialization points sampled straight from the priors before starting to '
                         'optimize the acquisition of points. Default = 262145 ')
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
                    default=300000,
                    help='Evidence points requested (including init-evidence). '
                         'Default = 300000 ')
    IO.add_argument('--data-dir',
                    help='Directory containing popPUNK distance files. ')
    IO.add_argument('--data-pref',
                    help='Prefix of popPUNK distance file(s). ')
    IO.add_argument('--seed',
                    type=int,
                    default=254,
                    help='Seed for random number generation. Default = 254. ')
    IO.add_argument('--summary',
                    choices=['quantile', 'mean'],
                    default="quantile",
                    help='Mode for summary statistics, either "mean" or "quantile". Default = "quantile". ')
    IO.add_argument('--mode',
                    choices=['ABC', 'BOLFI'],
                    default="BOLFI",
                    help='Mode for running model fit, either "ABC" or "BOLFI". Default = "BOLFI". ')
    IO.add_argument('--outpref',
                    default="./",
                    help='Output prefix. Default = "./"')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')

    return parser.parse_args()

class CustomPrior_s1(elfi.Distribution):
    def rvs(vec, size=1, random_state=None):
        locs = np.zeros(size)
        scales = 1 - vec
        t = scipy.stats.uniform.rvs(loc=locs, scale=scales, size=size, random_state=random_state)
        return t

class CustomPrior_s2(elfi.Distribution):
    def rvs(p1, p2, size=1, random_state=None):
        locs = np.zeros(size)
        scales = 1 - np.sum(np.array([p1, p2]), axis=0)
        t = scipy.stats.uniform.rvs(loc=locs, scale=scales, size=size, random_state=random_state)
        return t

class CustomPrior_s3(elfi.Distribution):
    def rvs(p1, p2, p3, size=1, random_state=None):
        locs = np.zeros(size)
        scales = 1 - np.sum(np.array([p1, p2, p3]), axis=0)
        t = scipy.stats.uniform.rvs(loc=locs, scale=scales, size=size, random_state=random_state)
        return t

class CustomPrior_s4(elfi.Distribution):
    def rvs(p1, p2, p3, p4, size=1, random_state=None):
        locs = np.zeros(size)
        scales = 1 - np.sum(np.array([p1, p2, p3, p4]), axis=0)
        t = scipy.stats.uniform.rvs(loc=locs, scale=scales, size=size, random_state=random_state)
        return t

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

def gen_distances_elfi(size_core, size_pan, prop_core_var, prop_acc_var, core_mu, acc_vs_core, avg_gene_freq, base_mu1, base_mu2,
                       base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4,
                       acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4,
                       batch_size=1, random_state=None):
    # determine vectors of core and accessory per-site mutation rates and variable regions
    core_mu_arr = np.array([core_mu] * batch_size)
    acc_vs_core = np.reshape(acc_vs_core, (-1, 1))
    acc_mu_arr = core_mu_arr * acc_vs_core
    core_var = np.round(size_core * prop_core_var)
    acc_var = np.round(size_pan * prop_acc_var)

    # ensure probabilities sum to 1, by taking sum and then determing overall proportion
    base_vec_sum = np.sum(np.array([base_mu1, base_mu2, base_mu3, base_mu4]), axis=0)
    core_vec_sum = np.sum(np.array([core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4]), axis=0)
    acc_vec_sum = np.sum(np.array([acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4]), axis=0)

    # generate vectors for mutation rates
    base_mu = np.stack((base_mu1, base_mu2, base_mu3, base_mu4), axis=1) / base_vec_sum[:,None]
    core_site = np.stack((core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4), axis=1) / core_vec_sum[:,None]
    acc_site = np.stack((acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4), axis=1) / acc_vec_sum[:,None]
    gene_mu = np.stack((1 - avg_gene_freq, avg_gene_freq), axis=1)

    core_site_mu = calc_man_vec(size_core, core_var, core_site)
    acc_site_mu = calc_man_vec(size_pan, acc_var, acc_site)

    # generate starting genomes
    core_ref = np.zeros((batch_size, size_core))
    acc_ref = np.zeros((batch_size, size_pan))
    for i in range(batch_size):
        core_ref[i] = np.random.choice([1, 2, 3, 4], size_core, p=base_mu[i])
        acc_ref[i] = np.random.choice([0, 1], size_pan, p=gene_mu[i])

    # mutate genomes
    core_query1 = sim_divergence_vec(core_ref, core_mu_arr, True, base_mu, core_site_mu)
    #core_query2 = sim_divergence_vec(core_ref, core_mu_arr, True, base_mu, core_site_mu)
    acc_query1 = sim_divergence_vec(acc_ref, acc_mu_arr, False, gene_mu, acc_site_mu)
    #acc_query2 = sim_divergence_vec(acc_ref, acc_mu_arr, False, gene_mu, acc_site_mu)

    # determine hamming and core distances
    hamming_core = np.zeros((batch_size, core_mu_arr.shape[1]))
    jaccard_acc = np.zeros((batch_size, core_mu_arr.shape[1]))

    for i in range(batch_size):
        for j in range(core_mu_arr.shape[1]):
            hamming_core[i][j] = distance.hamming(core_ref[i], core_query1[i][j])
            jaccard_acc[i][j] = distance.jaccard(acc_ref[i], acc_query1[i][j])

    # calculate euclidean distance to origin
    eucl = np.sqrt((hamming_core ** 2) + (jaccard_acc ** 2))

    return eucl

if __name__ == "__main__":
    #testing
    # size_core = 10000
    # size_pan = 1000
    # batch_size = 10000
    # N_samples = 100
    # qnt = 0.01
    # seed = 254
    # summary = "quantile"
    # data_dir = "/mnt/c/Users/sth19/PycharmProjects/PhD_project/distance_sim/distances"
    # data_pref = "GPSv4"
    # num_steps = 10
    # max_acc_vs_core = 1000
    # threads = 1
    # mode = "BOLFI"
    # outpref = "test_"
    # initial_evidence = 20
    # update_interval = 10
    # acq_noise_var = 0.1
    # n_evidence = 200
    # info_freq = 1000

    options = get_options()
    threads = options.threads
    data_dir = options.data_dir
    data_pref = options.data_pref
    size_core = options.core_size
    size_pan = options.pan_size
    batch_size = options.batch_size
    max_acc_vs_core = options.max_acc_vs_core
    num_steps = options.num_steps
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

    #set multiprocessing client
    elfi.set_client('multiprocessing')
    elfi.set_client(elfi.clients.multiprocessing.Client(num_processes=threads))

    # read in real files
    df = read_files(data_dir, data_pref)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(df["Core"].max())
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))

    # set constants
    # set evenly spaced core hamming values
    core_mu = np.linspace(0, max_real_core, num=num_steps)

    # set minimum prop_core_var and prop_acc_var based on number of sequence bins (hard coded at 5 at the moment)
    # min_prop_core_var = 5 / size_core
    # min_prop_acc_var = 5 / size_pan

    # set priors
    acc_vs_core = elfi.Prior('uniform', 0, max_acc_vs_core)
    avg_gene_freq = elfi.Prior('uniform', 0, 1)
    prop_core_var = elfi.Prior('uniform', 0, 1)
    prop_acc_var = elfi.Prior('uniform', 0, 1)

    # set priors based on remaining sum from previous allocations
    base_mu1 = elfi.Prior('uniform', 0, 1)
    base_mu2 = elfi.Prior('uniform', 0, 1)
    base_mu3 = elfi.Prior('uniform', 0, 1)
    base_mu4 = elfi.Prior('uniform', 0, 1)

    core_site_mu1 = elfi.Prior('uniform', 0, 1)
    core_site_mu2 = elfi.Prior('uniform', 0, 1)
    core_site_mu3 = elfi.Prior('uniform', 0, 1)
    core_site_mu4 = elfi.Prior('uniform', 0, 1)
    #core_site_mu5 = elfi.Prior('uniform', 0, 1)

    acc_site_mu1 = elfi.Prior('uniform', 0, 1)
    acc_site_mu2 = elfi.Prior('uniform', 0, 1)
    acc_site_mu3 = elfi.Prior('uniform', 0, 1)
    acc_site_mu4 = elfi.Prior('uniform', 0, 1)
    #acc_site_mu5 = elfi.Prior('uniform', 0, 1)

    #get observed data
    obs_core = df['Core'].to_numpy()
    obs_acc = df['Accessory'].to_numpy()

    # calculate euclidean distance to origin
    obs = np.sqrt((obs_core ** 2) + (obs_acc ** 2)).reshape(1, -1)

    Y = elfi.Simulator(gen_distances_elfi, size_core, size_pan, prop_core_var, prop_acc_var, core_mu, acc_vs_core,
                       avg_gene_freq, base_mu1, base_mu2, base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3,
                       core_site_mu4, acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4, observed=obs)

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

    if mode == "ABC":
        rej = elfi.Rejection(d, batch_size=batch_size, seed=seed)

        result = rej.sample(N_samples, quantile=qnt)
    else:
        log_d = elfi.Operation(np.log, d)
        bolfi = elfi.BOLFI(log_d, batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                           acq_noise_var=acq_noise_var, seed=seed)
        post = bolfi.fit(n_evidence=n_evidence)
        result = bolfi.sample(N_samples)

    with open(outpref + "ELFI_summary.txt", "w") as f:
        f.write(result)
