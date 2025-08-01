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
#from scipy.stats import wasserstein_distance_nd
import scipy.stats as ss
from KDE_distance import KDE_KL_divergence, KDE_JS_divergence, get_grid

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

def plot_negative_exponential(obs_df, outpref):
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

def safe_neg_inv(d):
    return -1 / np.log(np.clip(d, 1e-12, None))

def smooth_d(d):
    return np.power(np.clip(d, 0, None), 0.5)

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
    IO.add_argument('--epsilon',
                    type=float,
                    default=1e-5,
                    help='The minimum value for transformations to log-uniform space.  '
                         'Default = 1e-5 ')
    IO.add_argument('--noise-scale',
                    type=float,
                    default=1e-4,
                    help='The minimum value for transformations to log-uniform space.  '
                         'Default = 1e-4 ')
    IO.add_argument('--samples',
                    type=int,
                    default=100000,
                    help='No. samples for posterior estimation. Default = 100000 ')
    IO.add_argument('--init_evidence',
                    type=int,
                    default=50,
                    help='Number of initialization points sampled straight from the priors before starting to '
                         'optimize the acquisition of points. Default = 50 ')
    IO.add_argument('--threshold',
                    type=float,
                    default=None,
                    help='The threshold (bandwidth) for posterior  '
                         'Default = None ')
    IO.add_argument('--n_evidence',
                    type=int,
                    default=200,
                    help='Evidence points requested (including init-evidence). '
                         'Default = 200 ')
    IO.add_argument('--update-int',
                    type=int,
                    default=1,
                    help='Defines how often the GP hyperparameters are optimized. '
                            'Default = 1 ')
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
    
    # Add parameter specification arguments
    parser.add_argument('--param', nargs=4, action='append',
                        metavar=('NAME', 'MIN', 'MAX', 'DIST'),
                        help='Parameter to fit with BOLFI. Can be specified multiple times. DIST must be "uniform" or "loguniform"')
    
    # Fixed parameters (not fitted by BOLFI)
    parser.add_argument('--fixed-param', nargs=2, action='append',
                      metavar=('NAME', 'VALUE'),
                      help='Parameter with fixed value (not fitted by BOLFI). Can be specified multiple times.')

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
    """Prepare inputs for the simulator, handling both fixed and fitted parameters.
    
    This function processes the normalized parameter values [0,1] and transforms them
    back to their original ranges based on the parameter specifications.
    """

    # get sim_params
    fixed_params, other_params = inputs

    # Get the model from keyword arguments
    model = kwinputs.get('model')
    if model is None:
        raise ValueError("Model instance not found in keyword arguments")
    
    # Get parameter values in the order they were defined in the model
    param_values = {}
    param_index = 0
    
    # Get fixed parameters from keyword arguments
    #fixed_params = {k: v for k, v in kwinputs.items() if k not in ['model', 'meta', 'workdir', 'obs_file', 'epsilon']}
    
    # Process each node in the model (these are only the fitted parameters)
    for node_name, node in model.nodes.items():
        if isinstance(node, elfi.Prior) and node_name not in ['d', 'log_d', 'Y']:
            # Get the normalized parameter value [0,1]
            if param_index >= len(inputs):
                raise ValueError(f"Not enough input values for parameter {node_name}")
                
            norm_val = inputs[param_index][0] if isinstance(inputs[param_index], (list, np.ndarray)) else inputs[param_index]
            param_index += 1
            
            # Get parameter specification (default to uniform if not specified)
            spec = model.fitted_params[node_name]
            
            # Transform back to original range
            if spec['dist'] == 'uniform':
                # Linear scaling for uniform distribution
                real_val = spec['min'] + norm_val * (spec['max'] - spec['min'])
            else:  # loguniform
                # Log scaling for loguniform distribution
                real_val = from_unit_to_loguniform(norm_val, spec['min'], spec['max'])
            
            param_values[node_name] = real_val
    
    # Add fixed parameters to param_values
    param_values.update(fixed_params)
    param_values.update(other_params)

    # Add all parameters to kwinputs for the simulator
    for param_name, param_value in param_values.items():
        kwinputs[param_name] = param_value
    
    # Prepare a unique filename for parallel settings
    meta = kwinputs.get('meta', {})
    if 'workdir' in kwinputs and kwinputs['workdir'] is not None:
        filename = f"{kwinputs['workdir']}/{meta.get('model_name', 'model')}_{meta.get('batch_index', 0)}_{meta.get('submission_index', 0)}"
    else:
        filename = f"{meta.get('model_name', 'model')}_{meta.get('batch_index', 0)}_{meta.get('submission_index', 0)}"

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

    # read observations file
    obs = np.loadtxt(kwinputs['obs_file'], delimiter='\t', dtype='float64')

    divergence = KDE_JS_divergence(obs, simulations, eps=1e-12, log=True)

    # This will be passed to ELFI as the result of the command
    return divergence

if __name__ == "__main__":
    options = get_options()
    threads = options.threads
    obs_file = options.distfile
    N_samples = options.samples
    seed = options.seed
    outpref = options.outpref
    initial_evidence = options.init_evidence
    update_interval = options.update_int
    acq_noise_var = options.acq_noise_var
    n_evidence = options.n_evidence
    cluster = options.cluster
    load = options.load
    run_mode = options.run_mode
    max_distances = options.max_distances
    pansim_exe = options.pansim_exe
    chains = options.chains
    workdir = options.workdir
    threshold = options.threshold
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
    obs_df = np.loadtxt(obs_file, delimiter='\t', dtype='float64')

    # plot contours and curve fit
    plot_negative_exponential(obs_df, outpref)

    # set up model
    input_dim = 0
    m = elfi.ElfiModel(name='pansim_model')
    
    # Process fixed parameters first
    fixed_params = {}
    if hasattr(options, 'fixed_param') and options.fixed_param:
        for name, value in options.fixed_param:
            try:
                # remove quotation marks if string used
                value_parsed = value.replace('"', '')
                fixed_params[name] = float(value_parsed)
                print(f"Using fixed parameter: {name} = {str(fixed_params[name])}")
            except ValueError:
                raise ValueError(f"Invalid value for fixed parameter {name}: {value}. Must be a number.")
    
    # Check if we have at least one parameter to fit or fix
    if not hasattr(options, 'param') and not fixed_params:
        raise ValueError("No parameters specified. Use --param to specify parameters to fit "
                        "and/or --fixed-param to specify fixed parameters.")
    
    # Process user-specified parameters for BOLFI fitting
    fitted_params = {}
    if hasattr(options, 'param') and options.param:
        for name, min_val, max_val, dist in options.param:
            if name in fixed_params:
                raise ValueError(f"Parameter {name} cannot be both fixed (--fixed-param) and fitted (--param).")

            fitted_params[name] = {
                'min': float(min_val),
                'max': float(max_val),
                'dist': dist.lower()
            }
            
            # Validate distribution type
            if fitted_params[name]['dist'] not in ['uniform', 'loguniform']:
                raise ValueError(f"Invalid distribution type for parameter {name}. Must be 'uniform' or 'loguniform'.")
    
    # Add all user-specified parameters to the model (only those being fitted by BOLFI)
    m.fitted_params = {}
    for param_name, spec in fitted_params.items():
        # All parameters are defined as uniform in [0,1] space internally
        elfi.Prior('uniform', 0.0, 1.0, model=m, name=param_name)
        
        # Store parameter specifications as node data for later use
        m.fitted_params[param_name] = spec
        
        # Log the parameter being fitted
        print(f"Fitting parameter: {param_name} with {spec['dist']} prior in range [{spec['min']}, {spec['max']}]")

    # save observed parameters
    obs = np.array([0.0])

    # simulate and fit
    if run_mode == "sim":
        print("Simulating data...")

        # Define base command parameters with values from fixed_params and fitted_params
        base_params = []
        
        # Add fixed parameters to the command
        for param_name, param_value in fixed_params.items():
            base_params.append(f'--{param_name} {param_value}')
            
        # Add parameters being fitted to the command (they will be formatted later)
        for param_name in fitted_params.keys():
            base_params.append(f'--{param_name} {{{param_name}}}')
            
        # Add other required parameters that should always be passed
        other_params = {
            'seed': seed,
            'outpref': outpref,
            'max_distances': max_distances,
        }
        
        # Add other parameters that are not being fitted
        for param_name, param_value in other_params.items():
            if param_name not in fixed_params and param_name not in fitted_params:
                base_params.append(f'--{param_name} {param_value}')
            
        # Construct the final command
        command = pansim_exe + ' ' + ' '.join(base_params)

        WF_sim = elfi.tools.external_operation(command,
                                        prepare_inputs=prepare_inputs,
                                        process_result=process_result,
                                        stdout=False)

        WF_sim_vec = elfi.tools.vectorize(WF_sim)

        # save model
        save_path = outpref
        os.makedirs(save_path, exist_ok=True)      

        # Determine which parameters are being fitted
        bounds = {}
        
        # Get all parameters that are ELFI Priors (i.e., being fitted)
        for node_name, node in fitted_params.items():
            print(f"node_name, node: {node_name}, {node}")
            bounds[node_name] = (0.0, 1.0)
        
        print(bounds)

        # Create simulator with all required parameters
        other_other_params = {
            'workdir': workdir,
            'obs_file': obs_file,
            'epsilon': epsilon,
        }

        # add all non model parameters together
        other_params.update(other_other_params)
                
        # Create the simulator with all parameters
        elfi.Simulator(WF_sim_vec, fixed_params, other_params, name='sim', model=m, observed=obs)
        
        # Create array pool with all fitted parameters
        arraypool_fields = ['Y', 'd', 'log_d'] + [x for x in fitted_params.keys()]
        arraypool = elfi.ArrayPool(arraypool_fields, name="BOLFI_pool", prefix=save_path)

        m['sim'].uses_meta = True

        # Use euclidean distance between observed and simulated data
        elfi.Distance('euclidean', m['sim'], model=m, name='d')
        elfi.Operation(np.log, m['d'], model=m, name='log_d')
        
        # Set up the Gaussian Process model
        kernel = GPy.kern.Matern32(input_dim=len(bounds), ARD=True)
        # Alternative kernel: kernel = GPy.kern.RBF(input_dim=len(bounds), ARD=True)
        
        # Create target model with all parameters being fitted
        target_model = elfi.GPyRegression(
            parameter_names=list(bounds.keys()),
            bounds=bounds,
            kernel=kernel
        )
        
        # Initialize BOLFI
        mod = elfi.BOLFI(
            m['log_d'],
            batch_size=1,
            initial_evidence=initial_evidence,
            update_interval=update_interval,
            acq_noise_var=acq_noise_var,
            seed=seed,
            bounds=bounds,
            pool=arraypool,
            target_model=target_model
        )
        
        # Fit the model
        post = mod.fit(n_evidence=n_evidence, threshold=threshold)
        
        # Set up MCMC sampling
        sigma_proposals = {}
        for param_name, (low, high) in bounds.items():
            # Get parameter specification to determine if we should use log or linear scale for proposals
            node = m[param_name]
            
            # Scale proposals accordingly
            sigma_proposals[param_name] = (high - low) * covar_scaling
        
        # Run MCMC sampling
        result = mod.sample(
            N_samples,
            algorithm="metropolis",
            n_evidence=n_evidence,
            n_chains=chains,
            threshold=threshold,
            sigma_proposals=sigma_proposals
        )

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
        print("Loading models from: {}".format(load))
        if load is None:
            print('Error: Previously saved ELFI output required for "sample" mode. Please specify "--load".')
            sys.exit(1)

        try:
            # Parse the load path
            load_pref = os.path.dirname(load) if '/' in load else '.'
            load_name = os.path.basename(load)
            model_path = os.path.join(load_pref, "BOLFI_model")
            
            # Load the array pool
            print(f"Loading array pool from {load_pref}")
            arraypool = elfi.ArrayPool.open(name="BOLFI_pool", prefix=load_pref)
            print(f'Successfully loaded array pool with {len(arraypool)} batches')
            
            # Load the ELFI model
            print(f"Loading ELFI model from {model_path}")
            m = elfi.load_model(name="pansim_model", prefix=model_path)
            print("Successfully loaded ELFI model")
            
            # Get bounds from the model's parameter specifications
            bounds = {}
            for param_name, node in m.parameter_names.items():
                bounds[param_name] = (0.0, 1.0)
            
            if not bounds:
                raise ValueError("No parameter specifications found in the model. Cannot determine bounds.")
                
            print(f"Using parameter bounds from model: {bounds}")
            
        except Exception as e:
            print(f"Error loading model: {str(e)}")
            if 'arraypool' in locals():
                arraypool.close()
            raise
        
        kernel = GPy.kern.Matern32(input_dim=len(bounds), ARD=True)
        #kernel = GPy.kern.RBF(input_dim=len(bounds), ARD=True)
        target_model = elfi.GPyRegression(parameter_names=[x for x in bounds.keys()], bounds=bounds, kernel=kernel)
        mod = elfi.BOLFI(m['log_d'], batch_size=1, initial_evidence=initial_evidence, update_interval=update_interval,
                            acq_noise_var=acq_noise_var, seed=seed, bounds=bounds, pool=arraypool, target_model=target_model)

        post = mod.fit(n_evidence=n_evidence, threshold=threshold)
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