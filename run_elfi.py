import time

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import logging

import elfi
from simulate_divergence import *

def gen_distances_elfi(size_core, size_pan, prop_core_var, prop_acc_var, core_mu, acc_vs_core, avg_gene_freq, base_mu1, base_mu2,
                       base_mu3, core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4,
                       acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4,
                       batch_size=1, random_state=None):
    # ensure probabilities sum to 1
    vec_sum = np.sum(np.array([[base_mu1], [base_mu2], [base_mu3]]), axis=0)
    base_mu4 = 1 - vec_sum[0]
    vec_sum = np.sum(np.array([[core_site_mu1], [core_site_mu2], [core_site_mu3], [core_site_mu4]]), axis=0)
    core_site_mu5 = 1 - vec_sum[0]
    vec_sum = np.sum(np.array([[acc_site_mu1], [acc_site_mu2], [acc_site_mu3], [acc_site_mu4]]), axis=0)
    acc_site_mu5 = 1 - vec_sum[0]

    # generate vectors for mutation rates
    base_mu = np.concatenate((base_mu1, base_mu2, base_mu3, base_mu4), axis=1)
    core_site = np.concatenate((core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4, core_site_mu5), axis=1)
    acc_site = np.concatenate((acc_site_mu1, acc_site_mu2, acc_site_mu3, acc_site_mu4, acc_site_mu5), axis=1)
    gene_mu = np.concatenate((1 - avg_gene_freq, avg_gene_freq), axis=1, dtype=float)

    core_mu_arr = np.repeat(core_mu, repeats=batch_size, axis=0)
    acc_mu_arr = core_mu_arr * acc_vs_core

    core_var = np.round(size_core * prop_core_var)
    core_invar = size_core - core_var
    acc_var = np.round(size_pan * prop_acc_var)
    acc_invar = size_pan - acc_var

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
    core_query2 = sim_divergence_vec(core_ref, core_mu_arr, True, base_mu, core_site_mu)
    acc_query1 = sim_divergence_vec(acc_ref, acc_mu_arr, False, gene_mu, acc_site_mu)
    acc_query2 = sim_divergence_vec(acc_ref, acc_mu_arr, False, gene_mu, acc_site_mu)

    # determine hamming and core distances
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
    #test_vec = calc_man_vec(10000, np.array([[500], [1000]]), np.array([[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25]]))

    size_core = 10000
    size_pan = 1000
    prop_core_var = np.array([[0.25], [0.5]])
    prop_acc_var = np.array([[0.25], [0.5]])
    core_mu = np.array([[0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18]])
    acc_vs_core = np.array([[0.5], [1]])
    avg_gene_freq = np.array([[0.5], [0.75]])
    base_mu1 = np.array([[0.25], [0.7]])
    base_mu2 = np.array([[0.25], [0.1]])
    base_mu3 = np.array([[0.25], [0.1]])
    core_site_mu1 = np.array([[0.2], [0.6]])
    core_site_mu2 = np.array([[0.2], [0.1]])
    core_site_mu3 = np.array([[0.2], [0.1]])
    core_site_mu4 = np.array([[0.2], [0.1]])
    acc_site_mu1 = np.array([[0.2], [0.6]])
    acc_site_mu2 = np.array([[0.2], [0.1]])
    acc_site_mu3 = np.array([[0.2], [0.1]])
    acc_site_mu4 = np.array([[0.2], [0.1]])


    test_out = gen_distances_elfi(size_core, size_pan, prop_core_var, prop_acc_var, core_mu, acc_vs_core, avg_gene_freq, base_mu1,
                                  base_mu2, base_mu3, core_site_mu1, core_site_mu2, core_site_mu3, core_site_mu4, acc_site_mu1,
                                  acc_site_mu2, acc_site_mu3, acc_site_mu4, batch_size=2, random_state=None)


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