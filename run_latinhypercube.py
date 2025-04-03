import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
rng = np.random.default_rng()
import subprocess
import sys
from scipy.optimize import curve_fit
from scipy.stats import qmc

# run
"python run_latinhypercube.py --pan_genes 5000,50000 --core_genes 0,5000 --core_mu 0,0.1 --rate_genes1 0,100 --rate_genes2 0,100 --prop_genes2 0.0001,1.0 --prop_positive 0.001,1.0 --pos_lambda 0.00001,0.1 --pos_lambda 0.00001,0.1 --neg_lambda 0.00001,0.1 --avg_gene_freq 0.01,0.99 --HR_rate 0.0,10.0 --HGT_rate 0.0,10.0 --n_samples 10 --outpref sim_params"

def get_options():
    description = 'Run simulator of gene gain model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_elfi_simulator.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--pan_genes',
                    required=True,
                    help='Number of genes in pangenome, including core and accessory genes. ')
    IO.add_argument('--core_genes',
                    required=True,
                    help='Number of core genes in pangenome only.')
    IO.add_argument('--core_mu',
                    required=True,
                    help='Maximum pairwise distance for core genome.')
    IO.add_argument('--rate_genes1',
                    required=True,
                    help='Average number of accessory pangenome that mutates per generation in gene compartment 1. Must be >= 0.0.')
    IO.add_argument('--rate_genes2',
                    required=True,
                    help='Average number of accessory pangenome that mutates per generation in gene compartment 2. Must be >= 0.0.')
    IO.add_argument('--prop_genes2',
                    required=True,
                    help='Proportion of pangenome made up of compartment 2 genes. Must be 0.0 <= X <= 0.5.')
    IO.add_argument('--prop_positive',
                    required=True,
                    help='Proportion of pangenome made up of positively selected genes. Must be 0.0 <= X <= 1.0. If negative, neutral selection is simulated.')
    IO.add_argument('--pos_lambda',
                    required=True,
                    help='Lambda value for exponential distribution of positively selected genes. Must be > 0.0')
    IO.add_argument('--neg_lambda',
                    required=True,
                    help='Lambda value for exponential distribution of negatively selected genes. Must be > 0.0')
    IO.add_argument('--avg_gene_freq',
                    required=True,
                    help='Average gene frequency in accessory genome.')
    IO.add_argument('--HR_rate',
                    required=True,
                    help='Homologous recombination rate, as number of core sites transferred per core genome mutation.')
    IO.add_argument('--HGT_rate',
                    required=True,
                    help='HGT rate, as number of accessory sites transferred per core genome mutation.')
    IO.add_argument('--n_samples',
                    type=int,
                    default=100000,
                    help='Number of samples from hypercube.')
    IO.add_argument('--outpref',
                    default="sim",
                    help='Output prefix. Default = "sim"')
    # IO.add_argument('--seed',
    #                 type=int,
    #                 default=254,
    #                 help='Seed for random number generation. Default = 254. ')

    return parser.parse_args()

if __name__ == "__main__":
    options = get_options()
    core_mu = [float(x) for x in options.core_mu.split(",")]
    rate_genes1 = [float(x) for x in options.rate_genes1.split(",")]
    rate_genes2 = [float(x) for x in options.rate_genes2.split(",")]
    prop_genes2 = [float(x) for x in options.prop_genes2.split(",")]
    prop_positive = [float(x) for x in options.prop_positive.split(",")]
    pos_lambda = [float(x) for x in options.pos_lambda.split(",")]
    neg_lambda = [float(x) for x in options.neg_lambda.split(",")]
    #core_size = [float(x) for x in options.core_size.split(",")]
    pan_genes = [float(x) for x in options.pan_genes.split(",")]
    core_genes = [float(x) for x in options.core_genes.split(",")]
    avg_gene_freq = [float(x) for x in options.avg_gene_freq.split(",")]
    n_samples = options.n_samples
    outpref = options.outpref
    HR_rate = [float(x) for x in options.HR_rate.split(",")]
    HGT_rate = [float(x) for x in options.HGT_rate.split(",")]


    sampler = qmc.LatinHypercube(d=12, optimization="random-cd")

    # deal with negative numbers for prop_positive
    if prop_positive[0] != 0.0:
        prop_positive[0] *= -1

    l_bounds = [core_mu[0], rate_genes1[0], rate_genes2[0], prop_genes2[0], prop_positive[0], \
        pos_lambda[0], neg_lambda[0], pan_genes[0], core_genes[0], avg_gene_freq[0], HR_rate[0], HGT_rate[0]]
    u_bounds = [core_mu[1], rate_genes1[1], rate_genes2[1], prop_genes2[1], prop_positive[1], \
        pos_lambda[1], neg_lambda[1], pan_genes[1], core_genes[1], avg_gene_freq[1], HR_rate[1], HGT_rate[1]]
    
    #print(l_bounds)
    #print(u_bounds)
    
    sample = sampler.random(n=n_samples)
    scaled_sample = qmc.scale(sample, l_bounds, u_bounds)
    df = pd.DataFrame(scaled_sample)

    df.columns = ["core_mu", "rate_genes1", "rate_genes2", "prop_genes2", "prop_positive", \
        "pos_lambda", "neg_lambda", "pan_genes", "core_genes", "avg_gene_freq", "HR_rate", "HGT_rate"]
    
    df.pan_genes = df.pan_genes.round()
    df.core_genes = df.core_genes.round()

    df.to_csv(outpref + ".tsv", sep = "\t", index = False)

