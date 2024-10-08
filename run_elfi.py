import argparse

import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd
import elfi
import os
import sys
import matplotlib.pyplot as plt
import pickle
from scipy.spatial import distance
from scipy.optimize import curve_fit
from scipy.stats import wasserstein_distance_nd

# fit asymptotic curve using exponential decay
def asymptotic_curve(x, a, b, c):
    return a * (1 - np.exp(-b * x)) + c

# RMSE
def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))

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
    IO.add_argument('--pan_mu',
                    type=float,
                    default=1.0,
                    help='Diversification rate of accessory genome. Default = 1.0 (assumes complete saturation) ')
    IO.add_argument('--speed_fast',
                    type=float,
                    default=None,
                    help='Speed ratio at which a fast gene mutates over a slow gene. Default = None ')
    IO.add_argument('--n_gen',
                    type=int,
                    default=100,
                    help='Number of generations for Wright-Fisher model. Default = 100 ')
    IO.add_argument('--avg_gene_freq',
                    type=float,
                    default=0.5,
                    help='Average gene frequency in accessory genome.'
                         'Default = "0.5" ')
    IO.add_argument('--samples',
                    type=int,
                    default=1000,
                    help='No. samples for posterior estimation. Default = 1000 ')
    IO.add_argument('--init_evidence',
                    type=int,
                    default=5000,
                    help='Number of initialization points sampled straight from the priors before starting to '
                         'optimize the acquisition of points. Default = 5000 ')
    IO.add_argument('--n_evidence',
                    type=int,
                    default=5000,
                    help='Evidence points requested (including init-evidence). '
                         'Default = 5000 ')
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
    IO.add_argument('--chains',
                type=int,
                default=4,
                help='Number of chains for NUTS sampler. '
                        'Default = 4 ')
    IO.add_argument('--distfile',
                    required=True,
                    help='popPUNK distance file to fit to. ')
    IO.add_argument('--max_distances',
                    type=int,
                    default=100000,
                    help='Number of distances to sample with Pansim. Default = 100000')
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
            df = pd.read_csv(filename, index_col=None, header=None, sep=",")
        else:
            df = pd.read_csv(filename, index_col=None, header=None, sep="\t")

    if len(df.columns) == 2:
        df.rename(columns={df.columns[0]: "Core",
                           df.columns[1]: "Accessory"}, inplace=True)
    elif len(df.columns) == 4:
        # rename columns
        df.rename(columns={df.columns[0]: "Sample1", df.columns[1] : "Sample2", df.columns[2]: "Core",
                           df.columns[3]: "Accessory"}, inplace=True)
    else:
        print("Incorrect number of columns in distfile. Should be 2 or 4.")
        sys.exit(1)

    df['Core'] = pd.to_numeric(df['Core'])
    df['Accessory'] = pd.to_numeric(df['Accessory'])

    return df

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
    avg_gene_freq, pan_mu, proportion_fast, speed_fast, core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file = inputs
    
    # add to kwinputs
    kwinputs['avg_gene_freq'] = avg_gene_freq
    kwinputs['core_mu'] = core_mu
    kwinputs['pan_mu'] = pan_mu
    kwinputs['proportion_fast'] = proportion_fast
    kwinputs['speed_fast'] = speed_fast
    kwinputs['seed'] = seed
    kwinputs['pop_size'] = pop_size
    kwinputs['core_size'] = core_size
    kwinputs['pan_genes'] = pan_genes
    kwinputs['core_genes'] = core_genes
    kwinputs['n_gen'] = n_gen
    kwinputs['max_distances'] = max_distances
    kwinputs['obs_file'] = obs_file

    meta = kwinputs['meta']

    # Prepare a unique filename for parallel settings
    if workdir != None:
        filename = workdir + '/{model_name}_{batch_index}_{submission_index}'.format(**meta)
    else:
        filename = '{model_name}_{batch_index}_{submission_index}'.format(**meta)

    # Add the filenames to kwinputs
    kwinputs['output_filename'] = filename + '_out.txt'

    # Return new inputs that the command will receive
    return inputs, kwinputs

# Function to process the result of the simulation
def process_result(completed_process, *inputs, **kwinputs):
    output_filename = kwinputs['output_filename']

    # Read the simulations from the file.
    simulations = np.loadtxt(output_filename, delimiter='\t', dtype='float64')
    # Clean up the files after reading the data in
    os.remove(output_filename)

    obs_file = kwinputs['obs_file']
    obs = np.loadtxt(obs_file, delimiter='\t', dtype='float64')

    # get maximum core and accessory for distibution
    max_core = max(np.max(obs[:,0]), np.max(simulations[:,0]))
    max_acc = max(np.max(obs[:,1]), np.max(simulations[:,1]))
    #max_val = ma(max_core, max_acc)
    
    # process distributions
    sim_core = np.histogram(simulations[:,0], bins=200, range=(0, max_core))[0]
    sim_acc = np.histogram(simulations[:,1], bins=200, range=(0, max_acc))[0]
    sim_dist = (sim_core, sim_acc)

    obs_core = np.histogram(obs[:,0], bins=200, range=(0, max_core))[0]
    obs_acc = np.histogram(obs[:,1], bins=200, range=(0, max_acc))[0]
    obs_dist = (obs_core, obs_acc)

    js_core = js_distance(sim_dist, 0, obs_dist)
    js_pan = js_distance(sim_dist, 1, obs_dist)
    #was_dist = wasserstein_distance(sim_dist, obs_dist)

    dist = (js_core + js_pan) / 2
    #dist = was_dist
    #dist = js_pan

    # This will be passed to ELFI as the result of the command
    return dist

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
    pan_mu = options.pan_mu
    speed_fast = options.speed_fast
    workdir = options.workdir

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
    df = read_distfile(obs_file)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(df["Core"].max())
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))
    core_mu = max_real_core
    print("core_mu set to: {}".format(core_mu))

    # set up model
    m = elfi.ElfiModel(name='pansim_model')

    # set speed_fast
    if speed_fast == None or speed_fast < 1.0:
        max_value = 10 ** 6
        speed_fast = max_value
    
    print("speed_fast set to: {}".format(speed_fast))
    print("max pan_mu: {}".format(df["Accessory"].max()))

    # determine maximum divergence possible for pangenome
    pan_mu_upper = (pan_genes - core_genes) / pan_genes

    elfi.Prior('uniform', 0.0, (0.0 + pan_mu_upper), model=m, name='pan_mu')
    elfi.Prior('uniform', 0.0, 1.0, model=m, name='proportion_fast')

    #data = Y.generate(3)

    if run_mode == "sim":
        print("Simulating data...")
        bounds = {
            'pan_mu' : (0, pan_mu_upper),
            'proportion_fast' : (0, 1),
            #'speed_fast' : (0, max_value),
        }

        command = pansim_exe + ' --avg_gene_freq {avg_gene_freq} --pan_mu {pan_mu} --proportion_fast {proportion_fast} --speed_fast {speed_fast} --core_mu {core_mu} --seed {seed} --pop_size {pop_size} --core_size {core_size} --pan_genes {pan_genes} --core_genes {core_genes} --n_gen {n_gen} --max_distances {max_distances} --output {output_filename}'

        WF_sim = elfi.tools.external_operation(command,
                                        prepare_inputs=prepare_inputs,
                                        process_result=process_result,
                                        stdout=False)

        WF_sim_vec = elfi.tools.vectorize(WF_sim)

        elfi.Simulator(WF_sim_vec, avg_gene_freq, m['pan_mu'], m['proportion_fast'], speed_fast, core_mu, seed, pop_size, core_size, pan_genes, core_genes, n_gen, max_distances, workdir, obs_file, name='sim', model=m, observed=0)
        m['sim'].uses_meta = True

        # use cityblock as only single entry to calculate distance
        elfi.Distance('cityblock', m['sim'], model=m, name='d')
        elfi.Operation(np.log, m['d'], model=m, name='log_d')

        # save model
        save_path = outpref
        os.makedirs(save_path, exist_ok=True)
        arraypool = elfi.ArrayPool(['proportion_fast', 'Y', 'd', 'log_d', 'pan_mu'], name="BOLFI_pool", prefix=save_path)
        
        mod = elfi.BOLFI(m['log_d'], batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, bounds=bounds, pool=arraypool)

        #post = mod.fit(n_evidence=n_evidence)
        result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence, n_chains=chains)

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

        bounds = {
            'pan_mu' : (0, pan_mu_upper),
            'proportion_fast' : (0, 1),
            #'speed_fast' : (0, max_value),
        }
        mod = elfi.BOLFI(m['log_d'], batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, bounds=bounds, pool=arraypool)

        result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence, n_chains=chains)

        mod.plot_discrepancy()
        plt.savefig(outpref + "_BOLFI_discrepancy.png")
        plt.close()

        # # plot results
        # mod.plot_state()
        # plt.savefig(outpref + "_state.png")
        # plt.close()

        #plot MCMC traces
        result.plot_traces();
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

    sys.exit(0)