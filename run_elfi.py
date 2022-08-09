import time

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import logging

import elfi
from simulate_divergence import *

def gen_distances_elfi(size_core, size_pan, prop_core_var, prop_acc_var, core_mu, acc_vs_core, avg_gene_freq, base_mu1, base_mu2,
                       base_mu3, base_mu4, core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4, core_site_mu5,
                       acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4, acc_site_mu5,
                       batch_size=1, random_state=None):
    # generate vectors for mutation rates
    base_mu = np.transpose(np.array(base_mu1, base_mu2, base_mu3, base_mu4))
    core_site = np.transpose(np.array(core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4, core_site_mu5))
    acc_site = np.transpose(np.array(acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4, acc_site_mu5))

    core_mu_arr = np.repeat(core_mu, repeats=batch_size, axis=0)
    acc_mu_arr = core_mu_arr * acc_vs_core

    core_var = np.round(size_core * prop_core_var)
    core_invar = size_core - core_var
    acc_var = np.round(size_pan * prop_acc_var)
    acc_invar = size_pan - acc_var

    core_site_mu = calc_man_vec(core_var, core_site, size_core)
    acc_site_mu = calc_man_vec(acc_var, acc_site, size_pan)

    # mutate genomes
    core_query1, total_sites = sim_divergence(np.copy(core_var), core_mu, True, base_mu, core_site_mu)
    core_query2, total_sites = sim_divergence(np.copy(core_var), core_mu, True, base_mu, core_site_mu)
    acc_query1, total_sites = sim_divergence(np.copy(acc_var), acc_mu, False, avg_gene_freq, acc_site_mu)
    acc_query2, total_sites = sim_divergence(np.copy(acc_var), acc_mu, False, avg_gene_freq, acc_site_mu)

    # add core genes to accessory distances
    #acc_ref = np.append(acc_ref, np.ones(num_core))
    acc_query1 = np.append(acc_query1, np.ones(num_core))
    acc_query2 = np.append(acc_query2, np.ones(num_core))

    # add core invariant sites to core alignments
    #core_var = np.append(core_var, core_invar)
    core_query1 = np.append(core_query1, core_invar)
    core_query2 = np.append(core_query2, core_invar)

    acc_vs_core = 1
    pangenome_frac = 0

    if core_mu[index] != 0:
        # determine accessory vs core divergence rate
        prop_subs_core = (sites_mutated[0] + sites_mutated[1]) / (core_query1.size + core_query2.size)
        prop_subs_acc = (sites_mutated[2] + sites_mutated[3]) / (acc_query1.size + acc_query2.size)
        acc_vs_core = prop_subs_acc / prop_subs_core

    # determine pangenome_frac
    match = acc_query1 == acc_query2
    zeros_match = match[acc_query1 == 0]
    num_zero_match = np.count_nonzero(zeros_match)

    pangenome_frac = (acc_query1.size - num_zero_match) / acc_query1.size

    hamming_core = distance.hamming(core_query1, core_query2)
    hamming_acc = distance.hamming(acc_query1, acc_query2)
    jaccard_core = distance.jaccard(core_query1, core_query2)
    jaccard_acc = distance.jaccard(acc_query1, acc_query2)

    return (index, hamming_core, hamming_acc, jaccard_core, jaccard_acc, acc_vs_core, pangenome_frac)

def run_gen_distances():
    test = 1


def MA2(t1, t2, n_obs=100, batch_size=1, random_state=None):
    # Make inputs 2d arrays for numpy broadcasting with w
    t1 = np.asanyarray(t1).reshape((-1, 1))
    t2 = np.asanyarray(t2).reshape((-1, 1))
    random_state = random_state or np.random

    w = random_state.randn(batch_size, n_obs+2)  # i.i.d. sequence ~ N(0,1)
    x = w[:, 2:] + t1*w[:, 1:-1] + t2*w[:, :-2]
    return x

def autocov(x, lag=1):
    C = np.mean(x[:,lag:] * x[:,:-lag], axis=1)
    return C

if __name__ == "__main__":
    test_vec = calc_man_vec(10000, np.array([[500], [1000]]), np.array([[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25]]))


    # true parameters
    t1_true = 0.6
    t2_true = 0.2

    y_obs = MA2(t1_true, t2_true)

    # # Plot the observed sequence
    # plt.figure(figsize=(11, 6));
    # plt.plot(y_obs.ravel());
    #
    # # To illustrate the stochasticity, let's plot a couple of more observations with the same true parameters:
    # plt.plot(MA2(t1_true, t2_true).ravel());
    # plt.plot(MA2(t1_true, t2_true).ravel());

    # a node is defined by giving a distribution from scipy.stats together with any arguments (here 0 and 2)
    t1 = elfi.Prior(scipy.stats.uniform, 0, 2)

    # ELFI also supports giving the scipy.stats distributions as strings
    t2 = elfi.Prior('uniform', 0, 2)

    Y = elfi.Simulator(MA2, t1, t2, observed=y_obs)

    S1 = elfi.Summary(autocov, Y)
    S2 = elfi.Summary(autocov, Y, 2)  # the optional keyword lag is given the value 2

    # Finish the model with the final node that calculates the squared distance (S1_sim-S1_obs)**2 + (S2_sim-S2_obs)**2
    d = elfi.Distance('euclidean', S1, S2)

    test = 1