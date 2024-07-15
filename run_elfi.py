import argparse

import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd
import elfi
import os
import sys
import matplotlib.pyplot as plt

def get_options():
    description = 'Fit model to PopPUNK data using Approximate Baysesian computation'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_ELFI.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--run_mode',
                    required=True,
                    choices=['sim', 'sample'],
                    help='Which run mode to specify. Choices are "sim" or "sample".')
    IO.add_argument('--core_size',
                    type=int,
                    default=10000,
                    help='Number of positions in core genome. Default = 10000 ')
    IO.add_argument('--pan_size',
                    type=int,
                    default=10000,
                    help='Number of positions in pangenome. Default = 10000 ')
    IO.add_argument('--pop_size',
                    type=int,
                    default=1000,
                    help='Population size for Wright-Fisher model. Default = 1000 ')
    IO.add_argument('--n_gen',
                    type=int,
                    default=100,
                    help='Number of generations for Wright-Fisher model. Default = 100 ')
    IO.add_argument('--avg_gene_freq',
                    type=float,
                    default=0.5,
                    help='Average gene frequency in accessory genome. '
                         'Determines gene gain/loss rate e.g. 0.1 = gene gain/loss rate 1:9 '
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
    IO.add_argument('--distfile',
                    required=True,
                    help='popPUNK distance file to fit to. ')
    IO.add_argument('--max_distances',
                    type=int,
                    default=100000,
                    help='Number of distances to sample with Pansim. Default = 100000')
    IO.add_argument('--load',
                    default=None,
                    help='Directory of previous ELFI model. Required if running "sample" mode ')
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

# Function to prepare the inputs for the simulator. We will create filenames and write an input file.
def prepare_inputs(*inputs, **kwinputs):
    pan_mu, proportion_fast, speed_fast, core_mu, seed, pop_size, core_size, pan_size, n_gen, max_distances = inputs
    
    # add to kwinputs
    kwinputs['core_mu'] = core_mu
    kwinputs['pan_mu'] = pan_mu
    kwinputs['proportion_fast'] = proportion_fast
    kwinputs['speed_fast'] = speed_fast
    kwinputs['seed'] = seed
    kwinputs['pop_size'] = pop_size
    kwinputs['core_size'] = core_size
    kwinputs['pan_size'] = pan_size
    kwinputs['n_gen'] = n_gen
    kwinputs['max_distances'] = max_distances

    meta = kwinputs['meta']

    # Prepare a unique filename for parallel settings
    filename = '{model_name}_{batch_index}_{submission_index}'.format(**meta)

    # Add the filenames to kwinputs
    kwinputs['output_filename'] = filename + '_out.txt'

    # Return new inputs that the command will receive
    return inputs, kwinputs

# Function to process the result of the simulation
def process_result(completed_process, *inputs, **kwinputs):
    output_filename = kwinputs['output_filename']

    # Read the simulations from the file.
    simulations = np.loadtxt(output_filename, delimiter='\t', usecols=1, dtype='float64')
    obs_pan = np.histogram(simulations, bins=2000, range=(0, 1), density=True)[0]

    # Clean up the files after reading the data in
    os.remove(output_filename)

    # This will be passed to ELFI as the result of the command
    return obs_pan

if __name__ == "__main__":
    # #testing
    # core_size = 1200000
    # pan_size = 6000
    # #avg_gene_freq = 0.5
    # N_samples = 10
    # seed = 254
    # distfile = "GPSv4_distances_sample1.txt"
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
    distfile = options.distfile
    core_size = options.core_size
    pan_size = options.pan_size
    batch_size = options.batch_size
    N_samples = options.samples
    seed = options.seed
    outpref = options.outpref
    initial_evidence = options.init_evidence
    update_interval = options.update_int
    acq_noise_var = options.acq_noise_var
    n_evidence = options.n_evidence
    avg_gene_freq = options.avg_gene_freq
    cluster = options.cluster
    n_gen = options.ngen
    pop_size = options.pop_size
    load = options.load
    run_mode = options.run_mode
    max_distances = options.max_distances
    pansim_exe = options.pansim_exe

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
    df = read_distfile(distfile)

    # detemine highest core hamming distance, convert to real space using Jukes-Cantor
    max_hamming_core = float(df["Core"].max())
    max_jaccard_acc = float(df["Accessory"].max())
    max_real_core = (-3/4) * np.log(1 - (4/3 * max_hamming_core))
    core_mu = max_real_core

    #get observed data, normalise
    #obs_core = get_quantile(df['Core'].to_numpy())# / max_hamming_core)
    #obs_acc = get_quantile(df['Accessory'].to_numpy())# / max_jaccard_acc)
    obs_acc = np.histogram(df['Accessory'].to_numpy(), bins=2000, range=(0, 1), density=True)[0]

    # calculate euclidean distance to origin
    #obs = np.concatenate([obs_core, obs_acc])
    obs = obs_acc

    # set up model
    m = elfi.ElfiModel(name='pansim')

    # set priors
    max_value = 10 ** 6
    #elfi.Prior('uniform', 0, max_real_core, model=m, name='core_mu')
    elfi.Prior('uniform', 0.0, 1.0, model=m, name='pan_mu')
    elfi.Prior('uniform', 1.0, max_value, model=m, name='speed_fast')
    elfi.Prior('uniform', 0.0, 1.0, model=m, name='proportion_fast')

    command = pansim_exe + ' --pan_mu {pan_mu} --proportion_fast {proportion_fast} --speed_fast {speed_fast} --core_mu {core_mu} --seed {seed} --pop_size {pop_size} --core_size {core_size} --pan_size {pan_size} --n_gen {n_gen} --max_distances {max_distances} --output {output_filename}'

    WF_sim = elfi.tools.external_operation(command,
                                    prepare_inputs=prepare_inputs,
                                    process_result=process_result,
                                    stdout=False)
    
    WF_sim_vec = elfi.tools.vectorize(WF_sim)
    
    Y = elfi.Simulator(WF_sim_vec, m['pan_mu'], m['proportion_fast'], m['speed_fast'], core_mu, seed, pop_size, core_size, pan_size, n_gen, max_distances, observed=obs, name='sim')
    Y.uses_meta = True

    #data = Y.generate(3)

    d = elfi.Distance('jensenshannon', Y)
    log_d = elfi.Operation(np.log, d)

    if run_mode == "sim":
        print("Simulating data...")
        bounds = {
            'pan_mu' : (0, 1),
            'proportion_fast' : (0, 1),
            'speed_fast' : (0, max_value),
        }
        mod = elfi.BOLFI(log_d, batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, bounds=bounds)

        post = mod.fit(n_evidence=n_evidence)
        result = mod.sample(N_samples, algorithm="metropolis", info_freq=int(N_samples * 0.25))

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
        arraypool = elfi.OutputPool(['pan_mu', 'proportion_fast', 'speed_fast', 'Y', 'd'], prefix=save_path)
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

        log_d = elfi.Operation(np.log, d)
        bounds = {
            'pan_mu' : (0, 1),
            'proportion_fast' : (0, 1),
            'speed_fast' : (0, max_value),
        }
        mod = elfi.BOLFI(log_d, batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, bounds=bounds, pool=arraypool)

        result = mod.sample(N_samples, algorithm="metropolis", n_evidence=n_evidence, info_freq=int(N_samples * 0.25))

        mod.plot_discrepancy()
        plt.savefig(outpref + "_BOLFI_discrepancy.svg")
        plt.close()

        # plot results
        mod.plot_state()
        plt.savefig(outpref + "_state.svg")
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
    plt.savefig(outpref + '_marginals.svg')

    plt.clf
    plt.cla

    # plot paired marginals
    result.plot_pairs()
    plt.savefig(outpref + '_pairs.svg')
    plt.close()

    sys.exit(0)