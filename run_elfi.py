import argparse

import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd
import elfi
import GPy
import os
import sys
import matplotlib.pyplot as plt
import pickle
from scipy.spatial import distance
from scipy.optimize import curve_fit
from scipy.stats import wasserstein_distance_nd
import scipy.stats as ss

try:  # sklearn >= 0.22
    from sklearn.neighbors import KernelDensity
except ImportError:
    from sklearn.neighbors.kde import KernelDensity

# fit asymptotic curve using exponential decay
# b0 is asymptote, b1 is y-intercept, b2 is rate of decay
def negative_exponential(x, b0, b1, b2): # based on https://isem-cueb-ztian.github.io/Intro-Econometrics-2017/handouts/lecture_notes/lecture_10/lecture_10.pdf and https://www.statforbiology.com/articles/usefulequations/
    return b0 - (b0 - b1) * np.exp(-b2 * x)

def fit_negative_exponential(x, y, p0=[1.0, 0.0, 1.0], bounds=([0.0, 0.0, 0.0], [1.0, 1.0, np.inf])):
    return curve_fit(negative_exponential, x, y, p0=p0, bounds=bounds)

# Copyright John Lees and Nicholas Croucher 2025
def get_grid(minimum, maximum, resolution):
    x = np.linspace(minimum, maximum, resolution)
    y = np.linspace(minimum, maximum, resolution)
    xx, yy = np.meshgrid(x, y)
    xy = np.vstack([yy.ravel(), xx.ravel()]).T

    return(xx, yy, xy)

# Copyright John Lees and Nicholas Croucher 2025
def plot_scatter(X, out_prefix, x_fit, y_fit):
    # Plot results - max 1M for speed
    max_plot_samples = 1000000
    if X.shape[0] > max_plot_samples:
        X = utils.shuffle(X, random_state=random.randint(1,10000))[0:max_plot_samples,]

    # Kernel estimate uses scaled data 0-1 on each axis
    scale = np.amax(X, axis = 0)
    X /= scale

    plt.figure(figsize=(11, 8), dpi= 160, facecolor='w', edgecolor='k')
    xx, yy, xy = get_grid(0, 1, 100)

    # KDE estimate
    kde = KernelDensity(bandwidth=0.03, metric='euclidean',
                        kernel='epanechnikov', algorithm='ball_tree')
    kde.fit(X)
    z = np.exp(kde.score_samples(xy))
    z = z.reshape(xx.shape).T

    levels = np.linspace(z.min(), z.max(), 10)
    # Rescale contours
    plt.contour(xx*scale[0], yy*scale[1], z, levels=levels[1:], cmap='plasma')
    scatter_alpha = 1

    # Plot on correct scale
    plt.scatter(X[:,0]*scale[0].flat, X[:,1]*scale[1].flat, s=1, alpha=scatter_alpha)
    plt.plot(x_fit, y_fit, label=f"Negative exponential 3 param", color='red')

    plt.xlabel('Core distance (' + r'$\pi$' + ')')
    plt.ylabel('Accessory distance (' + r'$a$' + ')')
    plt.savefig(out_prefix + '_contours.png')
    plt.close()

# RMSE
def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))

# converts uniform[0,1] to logunifrom[min_val,max_val]
def from_unit_to_loguniform(u, min_val, max_val):
    return np.exp(np.log(min_val) + u * (np.log(max_val) - np.log(min_val)))

def to_normalised_log_uniform(x_real, low, high, eps=1e-12):
    x_safe = np.clip(x_real, low + eps, high)  # avoid log(0)
    log_low = np.log(low + eps)
    log_high = np.log(high)
    return (np.log(x_safe) - log_low) / (log_high - log_low)

def get_options():
    description = 'Fit model to PopPUNK data using Approximate Baysesian computation'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_ELFI.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--run_mode',
                    default='sim',
                    choices=['sim', 'sample'],
                    help='Which run mode to specify. Choices are "sim" or "sample".')
    IO.add_argument('--core_size',
                    type=int,
                    default=1200000,
                    help='Number of positions in core genome. Default = 1200000 ')
    IO.add_argument('--pan_genes',
                    type=int,
                    default=6000,
                    help='Number of genes in pangenome, including core and accessory genes. Default = 6000 ')
    IO.add_argument('--core_genes',
                    type=int,
                    default=2000,
                    help='Number of core genes in pangenome only. Default = 2000')
    IO.add_argument('--pop_size',
                    type=int,
                    default=1000,
                    help='Population size for Wright-Fisher model. Default = 1000 ')
    IO.add_argument('--n_gen',
                    type=int,
                    default=200,
                    help='Number of generations for Wright-Fisher model. Default = 200 ')
    IO.add_argument('--avg_gene_freq',
                    type=float,
                    default=0.5,
                    help='Average gene frequency in accessory genome.'
                         'Default = "0.5" ')
    IO.add_argument('--HR_rate',
                    type=float,
                    default=None,
                    help='Homologous recombination rate, as number of core sites transferred per core genome mutation.'
                         'If unspecified, will fit parameter. ')
    IO.add_argument('--HGT_rate',
                    type=float,
                    default=None,
                    help='HGT rate, as number of accessory sites transferred per core genome mutation.'
                         'If unspecified, will fit parameter. ')
    IO.add_argument('--recomb_max',
                    type=float,
                    default=10.0,
                    help='Maximum HGT and HR rate for parameterisation, as number of transfer events transferred per core genome mutation.'
                         'Default = 10.0. ')
    IO.add_argument('--competition',
                    action='store_true',
                    default=False,
                    help='Run simulator with competition.')
    IO.add_argument('--epsilon',
                    type=float,
                    default=1e-5,
                    help='The minimum value for transformations to log-uniform space.  '
                         'Default = 1e-5 ')
    IO.add_argument('--samples',
                    type=int,
                    default=100000,
                    help='No. samples for posterior estimation. Default = 100000 ')
    IO.add_argument('--init_evidence',
                    type=int,
                    default=1000,
                    help='Number of initialization points sampled straight from the priors before starting to '
                         'optimize the acquisition of points. Default = 1000 ')
    IO.add_argument('--threshold',
                    type=float,
                    default=None,
                    help='The threshold (bandwidth) for posterior  '
                         'Default = None ')
    IO.add_argument('--n_evidence',
                    type=int,
                    default=1000,
                    help='Evidence points requested (including init-evidence). '
                         'Default = 1000 ')
    IO.add_argument('--update-int',
                type=int,
                default=10,
                help='Defines how often the GP hyperparameters are optimized. '
                        'Default = 10 ')
    IO.add_argument('--acq-noise-var',
                type=float,
                default=0.01,
                help='Defines the diagonal covariance of noise added to the acquired points. '
                        'Default = 0.01 ')
    IO.add_argument('--chains',
                type=int,
                default=4,
                help='Number of chains for sampler. '
                        'Default = 4 ')
    IO.add_argument('--distfile',
                    required=True,
                    help='popPUNK distance file to fit to. ')
    IO.add_argument('--max_distances',
                    type=int,
                    default=100000,
                    help='Number of distances to sample with Pansim. Default = 100000')
    IO.add_argument('--covar-scaling',
                    type=float,
                    default=0.1,
                    help='Scaling of difference between lower and upper bounds of each parameter to be used for MCMC covariance. Default = 0.1')
    IO.add_argument('--load',
                    default=None,
                    help='Directory of previous ELFI model and pooled array, matching --outpref of previous run. Required if running "sample" mode ')
    IO.add_argument('--seed',
                    type=int,
                    default=254,
                    help='Seed for random number generation. Default = 254. ')
    IO.add_argument('--outpref',
                    default="PopPUNK-mod",
                    help='Output prefix. Default = "PopPUNK-mod"')
    IO.add_argument('--pansim_exe',
                    required=True,
                    help='Path to pansim executable.')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')
    IO.add_argument('--cluster',
                    action='store_true',
                    default=False,
                    help='Parallelise using ipyparallel if using cluster. Default = False')
    IO.add_argument('--workdir',
                default=None,
                help='Specify workdir to save intermediate files. If unset, will write to working directory.')
    

    return parser.parse_args()

def read_distfile(filename):
    # read first line, determine if csv
    with open(filename, "r") as f:
        first_line = f.readline()
        if "," in first_line:
            obs = pd.read_csv(filename, index_col=None, header=None, sep=",")
        else:
            obs = pd.read_csv(filename, index_col=None, header=None, sep="\t")

    if len(obs.columns) == 2:
        obs.rename(columns={obs.columns[0]: "Core",
                           obs.columns[1]: "Accessory"}, inplace=True)
    elif len(obs.columns) == 4:
        # rename columns
        obs.rename(columns={obs.columns[0]: "Sample1", obs.columns[1] : "Sample2", obs.columns[2]: "Core",
                           obs.columns[3]: "Accessory"}, inplace=True)
    else:
        print("Incorrect number of columns in distfile. Should be 2 or 4.")
        sys.exit(1)

    obs['Core'] = pd.to_numeric(obs['Core'])
    obs['Accessory'] = pd.to_numeric(obs['Accessory'])

    return obs

# process the summary statistic, 0 for core, 1 for accessory
def js_distance(sim, col, obs):
    y_sim = sim[col]
    y_obs = obs[col]

    #print("sim: {}".format(y_sim))
    #print("obs: {}".format(y_obs))

    js = distance.jensenshannon(y_obs, y_sim)
    #print("js: {}".format(js))
    
    return js

def wasserstein_distance(sim, obs):
    sim_array = np.column_stack(sim)
    obs_array = np.column_stack(obs)

    return wasserstein_distance_nd(obs_array, sim_array)

# Function to prepare the inputs for the simulator. We will create filenames and write an input file.
def prepare_inputs(*inputs, **kwinputs):
    avg_gene_freq, rate_genes1, rate_genes2, prop_genes2, core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file, HR_rate, HGT_rate, max_mu, recomb_max, epsilon = inputs
    
    # add to kwinputs
    kwinputs['avg_gene_freq'] = avg_gene_freq
    kwinputs['core_mu'] = core_mu
    kwinputs['rate_genes1'] = from_unit_to_loguniform(rate_genes1, epsilon, max_mu)
    # rate_genes2 is always faster make it so all genes in pangenome mutate once on average per generation
    kwinputs['rate_genes2'] = rate_genes2 * prop_genes2
    kwinputs['prop_genes2'] = from_unit_to_loguniform(prop_genes2, epsilon, 1.0)
    kwinputs['seed'] = seed
    kwinputs['pop_size'] = pop_size
    kwinputs['core_size'] = core_size
    kwinputs['pan_genes'] = pan_genes
    kwinputs['core_genes'] = core_genes
    kwinputs['n_gen'] = n_gen
    kwinputs['max_distances'] = max_distances
    kwinputs['HR_rate'] = from_unit_to_loguniform(HR_rate, epsilon, recomb_max)
    kwinputs['HGT_rate'] = from_unit_to_loguniform(HGT_rate, epsilon, recomb_max)

    meta = kwinputs['meta']

    # Prepare a unique filename for parallel settings
    if workdir != None:
        filename = workdir + '/{model_name}_{batch_index}_{submission_index}'.format(**meta)
        #filename = workdir + f'/{rate_genes1}_{prop_genes2}_{HGT_rate}'
    else:
        filename = '{model_name}_{batch_index}_{submission_index}'.format(**meta)
        #filename = f'/{rate_genes1}_{prop_genes2}_{HGT_rate}'

    # Add the filenames to kwinputs
    kwinputs['outpref'] = filename

    # Return new inputs that the command will receive
    return inputs, kwinputs

# Function to process the result of the simulation
def process_result(completed_process, *inputs, **kwinputs):
    output_filename = kwinputs['outpref'] + ".tsv"

    # Read the simulations from the file.
    simulations = np.loadtxt(output_filename, delimiter='\t', dtype='float64')
    # Clean up the files after reading the data in
    os.remove(output_filename)

    #obs_file = kwinputs['obs_file']
    #obs = np.loadtxt(obs_file, delimiter='\t', dtype='float64')

    # get maximum core and accessory for distibution
    #max_core = max(np.max(obs[:,0]), np.max(simulations[:,0]))
    #max_acc = max(np.max(obs[:,1]), np.max(simulations[:,1]))
    #max_val = ma(max_core, max_acc)
    
    # process distributions
    #sim_core = np.histogram(simulations[:,0], bins=500, range=(0, max_hamming_core))[0]
    #sim_acc = np.histogram(simulations[:,1], bins=500, range=(0, max_jacc_pan))[0]
    #sim = np.concatenate((sim_core, sim_acc), axis=0)

    try:
        popt, pcov = fit_negative_exponential(simulations[:,0], simulations[:,1])
        b0, b1, b2 = popt
        b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))
    except:
        b0, b1, b2, b0_err, b1_err, b2_err = -1e6, -1e6, -1e6, -1e6, -1e6, -1e6
    
    mean_acc = np.mean(simulations[:,1])

    #sim = np.array([b0, b1, b2, b0_err / (b0 + 1e-12), b1_err / (b1 + 1e-12), b2_err / (b2 + 1e-12)])
    #sim = np.array([b0, b1, b2])
    #sim = sim_acc
    #sim_dist = (sim_core, sim_acc)
    
    # based on random forest and hypercube (b0 for r1, mean_acc for prop2, b1 and b1_err for HGT)
    sim = np.array([b0, b2, b2_err, mean_acc])

    #obs_core = np.histogram(obs[:,0], bins=200, range=(0, max_core))[0]
    #obs_acc = np.histogram(obs[:,1], bins=200, range=(0, max_acc))[0]
    #obs_dist = (obs_core, obs_acc)

    #js_core = js_distance(sim_dist, 0, obs_dist)
    #js_pan = js_distance(sim_dist, 1, obs_dist)
    #was_dist = wasserstein_distance(sim_dist, obs_dist)

    #dist = (js_core + js_pan) / 2
    #dist = was_dist
    #dist = js_pan

    # This will be passed to ELFI as the result of the command
    return sim

if __name__ == "__main__":
    # #testing
    # core_size = 1200000
    # pan_genes = 6000
    # #avg_gene_freq = 0.5
    # N_samples = 10
    # seed = 254
    # obs_file = "GPSv4_distances_sample1.txt"
    # threads = 4
    # mode = "BOLFI"
    # outpref = "test"
    # initial_evidence = 20
    # update_interval = 10
    # acq_noise_var = 0.1
    # n_evidence = 200
    # info_freq = 1000
    # cluster = False
    # #complexity = "simple"
    # #schedule = "0.7,0.2,0.05"
    # pop_size = 1000
    # n_gen = 100
    # load = "test_pools/outputpool_254"
    # run_mode = "sim"
    # pansim_exe = "/Users/shorsfield/Documents/Software/Pansim/pansim/target/release/pansim"
    # max_distances = 100000

    options = get_options()
    threads = options.threads
    obs_file = options.distfile
    core_size = options.core_size
    pan_genes = options.pan_genes
    core_genes = options.core_genes
    N_samples = options.samples
    seed = options.seed
    outpref = options.outpref
    initial_evidence = options.init_evidence
    update_interval = options.update_int
    acq_noise_var = options.acq_noise_var
    n_evidence = options.n_evidence
    avg_gene_freq = options.avg_gene_freq
    cluster = options.cluster
    n_gen = options.n_gen
    pop_size = options.pop_size
    load = options.load
    run_mode = options.run_mode
    max_distances = options.max_distances
    pansim_exe = options.pansim_exe
    chains = options.chains
    workdir = options.workdir
    threshold = options.threshold
    HR_rate = options.HR_rate
    HGT_rate = options.HGT_rate
    competition = options.competition
    covar_scaling = options.covar_scaling
    epsilon = options.epsilon

    #set multiprocessing client
    os.environ['NUMEXPR_NUM_THREADS'] = str(threads)
    os.environ['NUMEXPR_MAX_THREADS'] = str(threads)

    if cluster == True:
        # must start ipyarallel cluster e.g. !ipcluster start -n threads --daemon
        elfi.set_client('ipyparallel')
    else:
        if threads > 1:
            elfi.set_client('multiprocessing')
            elfi.set_client(elfi.clients.multiprocessing.Client(num_processes=threads))
        else:
            elfi.set_client('native')

    # read in real files
    #obs = read_distfile(obs_file)
    obs_df = np.loadtxt(obs_file, delimiter='\t', dtype='float64')

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(np.max(obs_df[:,0]))
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))
    core_mu = max_real_core
    print("core_mu set to: {}".format(core_mu))

    # set up model
    input_dim = 0
    m = elfi.ElfiModel(name='pansim_model')

    # set max mutation rate to each gene being gained/lost once per generation, whole pangenome mutating for a single individual across the simulation
    max_mu = (pan_genes - core_genes) / n_gen
    elfi.Prior('uniform', 0.0, 1.0, model=m, name='rate_genes1')
    elfi.Prior('uniform', 0.0, 1.0, model=m, name='prop_genes2')
    
    # set rate_genes2 as total genome that can mutate, update with prop_genes2 to ensure each individual mutates 10x on average per generation (ensures saturation)
    rate_genes2 = (pan_genes - core_genes) * 10

    # set as arbitarily high value, 10 events per core genome mutation, or adjust to fit in normalised space
    recomb_max = options.recomb_max
    if HGT_rate == None:
        elfi.Prior('uniform', 0.0, 1.0, model=m, name='HGT_rate')
    else:
        HGT_rate = to_normalised_log_uniform(HGT_rate, 0.0, recomb_max, epsilon)
    
    if HR_rate == None:
        elfi.Prior('uniform', 0.0, 1.0, model=m, name='HR_rate')
    else:
        HR_rate = to_normalised_log_uniform(HR_rate, 0.0, recomb_max, epsilon)

    # fit negative_exponential curve
    popt, pcov = fit_negative_exponential(obs_df[:,0], obs_df[:,1])
    b0, b1, b2 = popt
    # get 1 std deviation error of parameters
    b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))

    mean_acc = np.mean(obs_df[:,1])

    # plot fit
    fig, ax = plt.subplots()
    ax.scatter(obs_df[:,0], obs_df[:,1], s=10, alpha=0.3)
    x_fit = np.linspace(0, obs_df[:,0].max(), 100)
    y_fit = negative_exponential(x_fit, *popt)
    ax.plot(x_fit, y_fit, label=f"Negative exponential 3 param", color='red')

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_annotate = 0.5 * xlim[1]  # 50% of the x-axis range
    y_annotate = 0.1 * ylim[1]  # 10% of the y-axis range

    # Calculate the initial rate at x=0
    print("Negative exponential 3 param, b0: {}, b1: {}, b2: {}".format(b0, b1, b2))

    ax.annotate("b0: {}, b1: {},\nb2: {}".format(round(b0, 3), round(b1, 3), round(b2, 3)), xy=(0, 0), xytext=(x_annotate, y_annotate),
             fontsize=10, color="green")

    ax.set_xlabel('Core distance (' + r'$\pi$' + ')')
    ax.set_ylabel('Accessory distance (' + r'$a$' + ')')

    fig.savefig(outpref + "_curve_fit.png")
    plt.close()

    plot_scatter(obs_df, outpref, x_fit, y_fit)

    # save observed parameters
    #obs = np.array([b0, b1, b2, b0_err / (b0 + 1e-12), b1_err / (b1 + 1e-12), b2_err / (b2 + 1e-12)])
    #obs = np.array([b0, b1, b2])
    obs = np.array([b0, b2, b2_err, mean_acc])

    # simulate and fit
    if run_mode == "sim":
        print("Simulating data...")

        if competition:
            command = pansim_exe + ' --avg_gene_freq {avg_gene_freq} --rate_genes1 {rate_genes1} --rate_genes2 {rate_genes2} --prop_genes2 {prop_genes2} --core_mu {core_mu} --seed {seed} --pop_size {pop_size} --core_size {core_size} --pan_genes {pan_genes} --core_genes {core_genes} --n_gen {n_gen} --max_distances {max_distances} --outpref {outpref} --HR_rate {HR_rate} --HGT_rate {HGT_rate} --competition'
        else:
            command = pansim_exe + ' --avg_gene_freq {avg_gene_freq} --rate_genes1 {rate_genes1} --rate_genes2 {rate_genes2} --prop_genes2 {prop_genes2} --core_mu {core_mu} --seed {seed} --pop_size {pop_size} --core_size {core_size} --pan_genes {pan_genes} --core_genes {core_genes} --n_gen {n_gen} --max_distances {max_distances} --outpref {outpref} --HR_rate {HR_rate} --HGT_rate {HGT_rate}'

        WF_sim = elfi.tools.external_operation(command,
                                        prepare_inputs=prepare_inputs,
                                        process_result=process_result,
                                        stdout=False)

        WF_sim_vec = elfi.tools.vectorize(WF_sim)

        # save model
        save_path = outpref
        os.makedirs(save_path, exist_ok=True)      

        if HR_rate != None and HGT_rate != None:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                #'rate_genes2' : (epsilon, max_mu),
                'prop_genes2' : (0.0, 1.0),
            }

            elfi.Simulator(WF_sim_vec, avg_gene_freq, m['rate_genes1'], rate_genes2, m['prop_genes2'], core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file, HR_rate, HGT_rate, max_mu, recomb_max, epsilon, name='sim', model=m, observed=obs)
            arraypool = elfi.ArrayPool(['Y', 'd', 'log_d', 'rate_genes1', 'prop_genes2'], name="BOLFI_pool", prefix=save_path)
        elif HR_rate != None and HGT_rate == None:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                'prop_genes2' : (0.0, 1.0),
                'HGT_rate' : (0.0, 1.0),
            }

            elfi.Simulator(WF_sim_vec, avg_gene_freq, m['rate_genes1'], rate_genes2, m['prop_genes2'], core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file, HR_rate, m['HGT_rate'], max_mu, recomb_max, epsilon, name='sim', model=m, observed=obs)
            arraypool = elfi.ArrayPool(['Y', 'd', 'log_d', 'rate_genes1', 'prop_genes2', 'HGT_rate'], name="BOLFI_pool", prefix=save_path)
        elif HR_rate == None and HGT_rate != None:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                'prop_genes2' : (0.0, 1.0),
                'HR_rate' : (0.0, 1.0),
            }
            elfi.Simulator(WF_sim_vec, avg_gene_freq, m['rate_genes1'], rate_genes2, m['prop_genes2'], core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file, m['HR_rate'], HGT_rate, max_mu, recomb_max, epsilon, name='sim', model=m, observed=obs)
            arraypool = elfi.ArrayPool(['Y', 'd', 'log_d', 'rate_genes1', 'prop_genes2','HR_rate'], name="BOLFI_pool", prefix=save_path)
        else:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                'prop_genes2' : (0.0, 1.0),
                'HR_rate' : (0.0, 1.0),
                'HGT_rate' : (0.0, 1.0),
            }
            elfi.Simulator(WF_sim_vec, avg_gene_freq, m['rate_genes1'], rate_genes2, m['prop_genes2'], core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file, m['HR_rate'], m['HGT_rate'], max_mu, recomb_max, epsilon, name='sim', model=m, observed=obs)
            arraypool = elfi.ArrayPool(['Y', 'd', 'log_d', 'rate_genes1', 'prop_genes2', 'HR_rate', 'HGT_rate'], name="BOLFI_pool", prefix=save_path)

        m['sim'].uses_meta = True

        # use euclidean between
        elfi.Distance('canberra', m['sim'], model=m, name='d')
        elfi.Operation(np.log, m['d'], model=m, name='log_d')     
        
        #kernel_ard = GPy.kern.RBF(input_dim=len(bounds), ARD=True, name='rbf')
        #target_model_ard = elfi.GPyRegression(parameter_names=[x for x in bounds.keys()], bounds=bounds, kernel=kernel_ard)
        mod = elfi.BOLFI(m['log_d'], batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, bounds=bounds, pool=arraypool)

        #post = mod.fit(n_evidence=n_evidence)
        result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence, n_chains=chains, threshold=threshold, sigma_proposals={key: (value[1] - value[0]) * covar_scaling for key, value in bounds.items()})

        mod.plot_discrepancy()
        plt.savefig(outpref + "_BOLFI_discrepancy.png")
        plt.close()

        # plot MCMC traces
        result.plot_traces(); 
        plt.savefig(outpref + '_BOLFI_traces.png')
        plt.close()

        with open(outpref + "_ELFI_summary.txt", "w") as f:
            print(result, file=f)
        
        arraypool.save()
        print('Files in', arraypool.path, 'are', os.listdir(arraypool.path))

        m.save(prefix=save_path + "/BOLFI_model")
        print('Model saved to ', save_path + "/BOLFI_model")

    else:
        print("Loading models in {}".format(load))
        if load == None:
            print('Previously saved ELFI output required for "sample" mode. Please specify "--load."')
            sys.exit(1)

        # parse filename
        load_pref = load.rsplit('/', 1)[0]
        load_name = load.rsplit('/', 1)[-1]

        arraypool = elfi.ArrayPool.open(name="BOLFI_pool", prefix=load_pref)
        print(arraypool[0])
        print('This pool has', len(arraypool), 'batches')

        m = elfi.load_model(name="pansim_model", prefix=load_pref + "/BOLFI_model")

        if HR_rate != None and HGT_rate != None:
                bounds = {
                    'rate_genes1' : (0.0, 1.0),
                    #'rate_genes2' : (epsilon, max_mu),
                    'prop_genes2' : (0.0, 1.0)
                }
        elif HR_rate != None and HGT_rate == None:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                #'rate_genes2' : (epsilon, max_mu),
                'prop_genes2' : (0.0, 1.0),
                'HGT_rate' : (0.0, recomb_max),
            }
        elif HR_rate == None and HGT_rate != None:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                #'rate_genes2' : (epsilon, max_mu),
                'prop_genes2' : (0.0, 1.0),
                'HR_rate' : (0.0, recomb_max),
            }
        else:
            bounds = {
                'rate_genes1' : (0.0, 1.0),
                #'rate_genes2' : (epsilon, max_mu),
                'prop_genes2' : (0.0, 1.0),
                'HR_rate' : (0.0, 1.0),
                'HGT_rate' : (0.0, 1.0),
            }
        
        #kernel_ard = GPy.kern.RBF(input_dim=len(bounds), ARD=True, name='rbf')
        #target_model_ard = elfi.GPyRegression(parameter_names=[x for x in bounds.keys()], bounds=bounds, kernel=kernel_ard)
        mod = elfi.BOLFI(m['log_d'], batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, pool=arraypool)

        result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence, n_chains=chains, threshold=threshold, sigma_proposals={key: (value[1] - value[0]) * covar_scaling for key, value in bounds.items()})

        mod.plot_discrepancy()
        plt.savefig(outpref + "_BOLFI_discrepancy.png")
        plt.close()

        # # plot results
        # mod.plot_state()
        # plt.savefig(outpref + "_state.png")
        # plt.close()

        #plot MCMC traces
        result.plot_traces()
        plt.savefig(outpref + '_BOLFI_traces.png')
        plt.close()

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
    plt.close()

    # generate output
    # Define real-world log-uniform min/max for each parameter
    param_names = mod.model.parameter_names
    param_bounds = {
        'rate_genes1': (epsilon, max_mu),
        'prop_genes2': (epsilon, 1.0),
        'HR_rate': (epsilon, recomb_max),
        'HGT_rate': (epsilon, recomb_max),
        # Add all your parameter bounds here
    }

    # Extract normalized samples
    X = mod.target_model.X  # shape: (n_samples, n_params)
    Y = mod.target_model.Y.flatten()

    # Transform normalized samples back to real-world scale
    X_real = np.zeros_like(X)
    df_post = pd.DataFrame(result.samples)

    summary_rows = []

    for i, pname in enumerate(param_names):
        if pname in param_bounds:
            min_val, max_val = param_bounds[pname]
            X_real[:, i] = from_unit_to_loguniform(X[:, i], min_val, max_val)
            df_post[pname] = from_unit_to_loguniform(df_post[pname], min_val, max_val)

            # Compute posterior summaries
            mean = np.mean(df_post[pname])
            ci_2_5 = np.percentile(df_post[pname], 2.5)
            ci_97_5 = np.percentile(df_post[pname], 97.5)
            median = np.median(df_post[pname])

            summary_rows.append({
                'parameter': pname,
                'mean': mean,
                'median': median,
                'CI_2.5%': ci_2_5,
                'CI_97.5%': ci_97_5
            })

    # Store in DataFrame and save
    df_evidence_scaled = pd.DataFrame(X_real, columns=param_names)
    df_evidence_scaled['discrepancy'] = Y
    df_evidence_scaled.to_csv(outpref + '_gp_evidence.csv', index=False)        
    df_post.to_csv(outpref + '_mcmc_posterior_samples.csv', index=False)
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(outpref + '_parameter_estimates_summary.csv', index=False)

    sys.exit(0)