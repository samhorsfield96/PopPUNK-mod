import math

import pandas as pd
from scipy.spatial import distance
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import random
from math import e
from scipy.stats import gamma
from numba import jit

def recurse_prob(x, weight):
    if x == 1:
        return weight
    else:
        return (1 - recurse_prob(x - 1, weight)) * weight

def calc_man_vec(array_size, vec_size, bin_probs):
    no_split = np.shape(bin_probs)[1]

    # get bins for pdf
    integer = vec_size // no_split
    mod = np.reshape((vec_size % no_split).astype(np.int64), (-1, 1))

    bins = np.transpose(np.array([integer] * np.shape(bin_probs)[1]).astype(int))

    # add modulus to bins to assign all sites
    bins[:, 0] += mod[:, 0]

    # assign probabilities to sites
    site_mu = np.zeros((np.shape(bin_probs)[0], array_size))

    scaled_bin_probs = bin_probs / bins

    start = np.zeros(np.shape(bin_probs)[0], dtype=np.int64)
    end = np.copy(bins[:, 0])

    r = np.arange(site_mu.shape[1])

    # for testing
    #non_zero = np.zeros(np.shape(bin_probs)[0])

    for index in range(no_split):
        mask = (start[:, None] <= r) & (end[:, None] > r)

        scaled_bins_index = scaled_bin_probs[:, index]

        for entry in range(scaled_bins_index.size):
            site_mu[entry][mask[entry]] = scaled_bins_index[entry]

        if index < np.shape(bin_probs)[1] - 1:
            start += bins[:, index]
            end += bins[:, index + 1]

        # for testing
        #non_zero += np.count_nonzero(mask, axis=1)

    #for testing
    # for row in range(site_mu.shape[0]):
    #     sum_sites_mu = np.sum(site_mu[row])

    return site_mu

@jit(nopython=True)
def sim_divergence_vec(ref, mu, core, freq, site_mu):
    num_sites_vec = np.empty_like(mu)
    np.round(ref.shape[1] * mu, 0, num_sites_vec)
    num_sites_vec = num_sites_vec.astype(np.int64)

    if core:
        choices = np.array([1, 2, 3, 4])
    else:
        choices = np.array([0, 1])

    index_array = np.arange(ref[0].size)

    #total_sites_vec = np.zeros((mu.shape[0]))

    # create 3d array, columns are each base, rows are each mu rate, depth is num batches
    query = np.zeros((mu.shape[0], mu.shape[1], ref.shape[1]))

    # iterate until all required sites mutated for given mutation rate
    for i in range(mu.shape[0]):
        for j in range(mu.shape[1]):
            query[i][j] = ref[i].copy()

            num_sites = num_sites_vec[i][j]
            if num_sites > 0:
                total_sites = 0
                while total_sites < num_sites:
                    to_sample = num_sites - total_sites

                    # pick all sites to be mutated
                    sites = index_array[np.searchsorted(np.cumsum(site_mu[i]), np.random.rand(to_sample))]

                    # determine number of times each site can be mutated
                    sample_sites = np.array(list(set(sites)))

                    # determine sites with and without change
                    changes = choices[np.searchsorted(np.cumsum(freq[i]), np.random.rand(sample_sites.size))]

                    non_mutated = query[i][j][sample_sites] == changes

                    query[i][j][sample_sites] = changes

                    # determine number of actual changes
                    total_sites += changes.size - np.count_nonzero(non_mutated)

    return query


def calc_man(vec_size, bin_probs):
    no_split = len(bin_probs)

    # get bins for pdf
    i, d = divmod(vec_size, no_split)
    mod = vec_size % no_split

    bins = [i for x in range(no_split)]

    # add modulus to bins to assign all sites
    for i in range(mod):
        bins[i] += 1

    # assign probabilities to sites
    site_mu = np.zeros(vec_size)

    start = 0
    for index, prob in enumerate(bin_probs):
        scaled_prob = prob / bins[index]
        site_mu[start : start + bins[index]] = scaled_prob

        start += bins[index]

    #sum_sites_mu = np.sum(site_mu)

    return site_mu

def calc_gamma(vec_size, no_split, shape, scale, sim_index):
    split_range = np.linspace(gamma.ppf(0.001, shape),
                    gamma.ppf(0.999, shape), no_split)
    cdf = gamma.cdf(split_range, shape, scale)

    # determine total culumative probability of bin
    bin_probs = [None] * no_split
    bin_probs[0] = cdf[0]
    for i in range(1, cdf.size - 1):
        bin_probs[i] = cdf[i] - cdf[i - 1]

    # ensure all add up to 1
    if no_split > 1:
        bin_probs[-1] = 1 - cdf[-2]
    else:
        bin_probs[0] = 1.0

    if sim_index == 0:
        print(bin_probs)

    #sum_bins_mu = np.sum(bin_probs)

    # get bins for pdf
    i, d = divmod(vec_size, no_split)
    mod = vec_size % no_split

    bins = [i for x in range(no_split)]

    # add modulus to bins to assign all sites
    for i in range(mod):
        bins[i] += 1

    # assign probabilities to sites
    site_mu = np.zeros(vec_size)

    start = 0
    for index, prob in enumerate(bin_probs):
        scaled_prob = prob / bins[index]
        site_mu[start : start + bins[index]] = scaled_prob

        start += bins[index]

    #sum_sites_mu = np.sum(site_mu)

    return site_mu

def sim_divergence(ref, mu, core, freq, site_mu):
    num_sites = round(len(ref) * mu)

    if core:
        choices = np.array([1, 2, 3, 4])
    else:
        choices = np.array([0, 1])
        freq = [1 - freq, freq]

    index_array = np.arange(ref.size)

    # create 3d array, columns are each base, rows are each mu rate, depth is num batches
    query = ref.copy()
    total_sites = 0
    if num_sites > 0:
        while total_sites < num_sites:
            to_sample = num_sites - total_sites

            # pick all sites to be mutated
            sites = np.random.choice(index_array, to_sample, p=site_mu)

            # determine number of times each site can be mutated
            sample_sites = np.array(list(set(sites)))

            # determine sites with and without change
            changes = np.random.choice(choices, sample_sites.size, p=freq)

            non_mutated = query[sample_sites] == changes

            query[sample_sites] = changes

            # determine number of actual changes
            total_sites += changes.size - np.count_nonzero(non_mutated)

    return query, total_sites

def gen_distances(index, core_var, acc_var, core_invar, num_core, core_mu, acc_mu, adj, avg_gene_freq, base_mu,
                  core_site_mu, acc_site_mu):
    sites_mutated = []

    # mutate genomes
    core_query1, total_sites = sim_divergence(core_var, core_mu[index], True, base_mu, core_site_mu)
    sites_mutated.append(total_sites)
    core_query2, total_sites = sim_divergence(core_var, core_mu[index], True, base_mu, core_site_mu)
    sites_mutated.append(total_sites)
    acc_query1, total_sites = sim_divergence(acc_var, acc_mu[index], False, avg_gene_freq, acc_site_mu)
    sites_mutated.append(total_sites)
    acc_query2, total_sites = sim_divergence(acc_var, acc_mu[index], False, avg_gene_freq, acc_site_mu)
    sites_mutated.append(total_sites)

    if adj == True:
        # add core genes to accessory distances
        #acc_ref = np.append(acc_ref, np.ones(num_core))
        acc_query1 = np.append(acc_query1, np.ones(num_core))
        acc_query2 = np.append(acc_query2, np.ones(num_core))

        # add core invariant sites to core alignments
        #core_var = np.append(core_var, core_invar)
        core_query1 = np.append(core_query1, core_invar)
        core_query2 = np.append(core_query2, core_invar)

    acc_vs_core = 1

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

def model(x, c0, c1, c2, c3):
    return (1/2 * (1 - np.sqrt((1 - (4/3 * x)) ** (3 * c0)))) / (c1 + c2 * x + c3 * (x ** 2))

def model2(x, c0, c1, c2):
    #model 1 works
    #return c0 + c1 * x - c2 * np.exp(-c3 * x)
    # model 2 doesn't work
    #c0 - c1 * np.exp(-c2 * x)
    # model 3 doesn't work
    #c0 - c1 * c2 * np.exp(-c3 * x)
    # model 4
    return c0 + c1 * x + c2 * (x ** 2)

def fit_cvsa_curve(hamming_core_sim, jaccard_accessory_sim):
    sim = 0
    reg_x = np.zeros(len(hamming_core_sim) * len(hamming_core_sim[0]))
    reg_y = np.zeros(len(jaccard_accessory_sim) * len(jaccard_accessory_sim[0]))

    for hamming_core, jaccard_accessory in zip(hamming_core_sim, jaccard_accessory_sim):
        # make data, as comparing two diverged sequences, multiply by 2
        x = np.array(hamming_core)
        y = np.array(jaccard_accessory)

        reg_x[sim * x.size : (sim * x.size) + x.size] = x
        reg_y[sim * y.size : (sim * y.size) + y.size] = y

        sim += 1

    #reg_x = reg_x.reshape((-1, 1))
    try:
        c, cov = curve_fit(model, reg_x, reg_y, maxfev=5000, bounds=(([0,0,0,-1]), (np.inf,1,1,0)))
    except RuntimeError:
        c = np.zeros(4)

    return c

def check_panfrac(distances, pangenome_fracs, outpref):
    reg_x = np.zeros(len(distances[0]) * len(distances[0][0]))
    reg_y = np.zeros(len(distances[1]) * len(distances[1][0]))

    for var1, var2, name in zip(distances, pangenome_fracs, ("hamming_core", "jaccard_accessory")):
        # plot
        fig, ax = plt.subplots()

        for j in range(len(var1)):
            x = np.array(var1[j])
            y = np.array(var2[j])

            reg_x[(j) * x.size: ((j) * x.size) + x.size] = x
            reg_y[(j) * y.size: ((j) * y.size) + y.size] = y

            ax.plot(x, y, linewidth=2.0, label="Sim" + str(j + 1))

        if name == "hamming_core":
            try:
                c, cov = curve_fit(model2, reg_x, reg_y, maxfev=5000, bounds=(([0,0,-np.inf]), (1,np.inf,0)))

                with open(outpref + "pangenome_frac_model_parameters.txt", "w") as f:
                    f.write(np.array2string(c))

                x = np.array(distances[0][0])
                y = np.array([model2(j, c[0], c[1], c[2]) for j in x])
                ax.plot(x, y, linewidth=2.0, label="Model")
                print("Model parameters for hamming core vs. pangenome frac:")
                print(c)
            except RuntimeError:
                pass

        ax.set_xlabel(name)
        ax.set_ylabel("pangenome fraction")
        ax.legend()

        fig.savefig(outpref + name + "_vs_pangenome_frac.png")

        ax.clear()
        plt.clf
        plt.cla

def generate_graph(mu_rates, distances, mu_names, distance_names, outpref, core_adj, acc_adj, adjusted):
    for var1, var2, name1, name2 in zip(mu_rates, distances, mu_names, distance_names):
        # plot
        fig, ax = plt.subplots()

        # make data, as comparing two diverged sequences, multiply by 2
        x = np.array(var1) * 2
        for j in range(len(var2)):
            y = np.array(var2[j])

            ax.plot(x, y, linewidth=2.0, label="Sim" + str(j + 1))

        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]

        # plot Jukes-cantor relationship
        if "Core" in name1:
            if adjusted:
                y = (3/4 * (1 - e ** (-(4/3) * x))) / core_adj
            else:
                y = 3 / 4 * (1 - e ** (-(4 / 3) * x))
        else:
            if adjusted:
                y = (1/2 * (1 - e ** (-(2/1) * x))) / acc_adj
            else:
                y = 1 / 2 * (1 - e ** (-(2 / 1) * x))
        ax.plot(x, y, linewidth=2.0, label="Jukes-Cantor")

        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0, label="y=x")
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(0, 1)
        ax.set_xlabel(name1)
        ax.set_ylabel(name2)
        ax.legend()

        fig.savefig(outpref + name2 + ".png")

        ax.clear()
        plt.clf
        plt.cla

    # plot core hamming vs. accessory jaccard, predict relationship
    fig, ax = plt.subplots()
    sim = 1
    reg_x = np.zeros(len(distances[0]) * len(distances[0][0]))
    reg_y = np.zeros(len(distances[3]) * len(distances[3][0]))
    for hamming_core, jaccard_accessory in zip(distances[0], distances[3]):
        # make data, as comparing two diverged sequences, multiply by 2
        x = np.array(hamming_core)
        y = np.array(jaccard_accessory)

        ax.plot(x, y, linewidth=2.0, label="Sim" + str(sim))

        reg_x[(sim - 1) * x.size : ((sim - 1) * x.size) + x.size] = x
        reg_y[(sim - 1) * y.size : ((sim - 1) * y.size) + y.size] = y

        sim += 1

    #fit model, determine uncertainty
    c, cov = curve_fit(model, reg_x, reg_y, maxfev=5000, bounds=(([0,0,0,-np.inf]), (np.inf,1,np.inf,0)))
    d_c0 = np.sqrt(cov[0][0])
    d_c1 = np.sqrt(cov[1][1])
    d_c2 = np.sqrt(cov[2][2])
    d_c3 = np.sqrt(cov[3][3])

    # predict using new model
    x = np.array(distances[0][0])
    y = np.array([model(j, c[0], c[1], c[2], c[3]) for j in x])
    ax.plot(x, y, linewidth=2.0, label="Model")
    print("Accessory vs. core model parameters:")
    print(c)
    print("[" + str(d_c0) + " " + str(d_c1) + " " + str(d_c2) + " " + str(d_c3) + "]")

    with open(outpref + "model_parameters.txt", "w") as f:
        f.write(np.array2string(c))
        f.write("[" + str(d_c0) + " " + str(d_c1) + " " + str(d_c2) + " " + str(d_c3) + "]")

    lims = [
        np.min([ax.get_xlim(), ax.get_xlim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_xlim()]),  # max of both axes
    ]

    # ax.set_ylim(0, 1)
    # ax.set_xlim(lims)
    ax.set_xlabel("Core Hamming distance")
    ax.set_ylabel("Accessory Jaccard Distance")
    ax.legend()

    fig.savefig(outpref + "core_vs_acc.png")
    plt.cla

    plt.clf

    # plot accessory jaccard vs. accessory hamming
    fig, ax = plt.subplots()
    sim = 1
    reg_x = np.zeros(len(distances[1]) * len(distances[1][0]))
    reg_y = np.zeros(len(distances[3]) * len(distances[3][0]))
    for hamming, jaccard in zip(distances[1], distances[3]):
        # make data, as comparing two diverged sequences, multiply by 2
        x = np.array(hamming)
        y = np.array(jaccard)
        reg_x[(sim - 1) * x.size : ((sim - 1) * x.size) + x.size] = x
        reg_y[(sim - 1) * y.size : ((sim - 1) * y.size) + y.size] = y

        ax.plot(x, y, linewidth=2.0, label="Sim " + str(sim))
        sim += 1

    # try:
    #     c, cov = curve_fit(model2, reg_x, reg_y, maxfev=5000)
    #
    #     with open(outpref + "acc_hamming_vs_jacc_parameters.txt", "w") as f:
    #         f.write(np.array2string(c))
    #
    #     x = np.array(distances[1][0])
    #     y = np.array([model2(j, c[0], c[1]) for j in x])
    #     ax.plot(x, y, linewidth=2.0, label="Model")
    #     print("Model parameters for hamming acc vs. jaccard acc:")
    #     print(c)
    # except RuntimeError:
    #     pass

    ax.set_xlabel("Accessory Hamming distance")
    ax.set_ylabel("Accessory Jaccard Distance")
    ax.legend()

    fig.savefig(outpref + "acc_hamming_vs_jaccard.png")
    plt.cla

    plt.clf

#for testing
if __name__ == "__main__":
    threads = 4
    size_core = 1140000
    num_core = 1194
    num_pangenome = 5442
    num_sim = 2
    core_num_var = 106196
    core_sites = 5
    acc_sites = 5

    # core_sites_man = [0.7, 0.1, 0.1, 0.1]
    # acc_sites_man = [0.7, 0.1, 0.1, 0.1]

    core_sites_man = None
    acc_sites_man = None

    # base frequencies are alphabetical, A, C, G, T
    base_freq = [0.25, 0.25, 0.25, 0.25]
    base_choices = [1, 2, 3, 4]

    base_mu = [0.25, 0.25, 0.25, 0.25]

    # gene presence/absence frequencies are in order 0, 1
    gene_freq = [0.5, 0.5]
    gene_choices = [0, 1]
    avg_gene_freq = 0.5

    # core mu is number of differences per base of alignment
    core_mu = [0.1 * i for i in range(0, 21, 2)]

    # acc mu is number of differences per gene in accessory genome
    acc_mu = [0.1 * i for i in range(0, 11)]

    adj = True

    hamming_core = [None] * len(core_mu)
    hamming_acc = [None] * len(core_mu)
    jaccard_core = [None] * len(core_mu)
    jaccard_acc = [None] * len(core_mu)

    # generate references
    core_ref = np.random.choice(base_choices, size_core, p=base_freq)
    size_acc = num_pangenome - num_core
    acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)

    hamming_core_sims = [None] * num_sim
    hamming_acc_sims = [None] * num_sim
    jaccard_core_sims = [None] * num_sim
    jaccard_acc_sims = [None] * num_sim

    # with Pool(processes=threads) as pool:
    #     for i in range(num_sim):
    #         print("Simulation " + str(i + 1))
    #
    #         hamming_core = [None] * len(core_mu)
    #         hamming_acc = [None] * len(acc_mu)
    #         jaccard_core = [None] * len(core_mu)
    #         jaccard_acc = [None] * len(acc_mu)
    #
    #
    #         for ind, hcore, hacc, jcore, jacc in tqdm.tqdm(pool.imap(
    #                 partial(gen_distances, core_ref=core_ref, acc_ref=acc_ref, core_invar=core_invar,
    #                         num_core=num_core, core_mu=core_mu, acc_mu=acc_mu),
    #                 range(0, len(core_mu))), total=len(core_mu)):
    #             hamming_core[ind] = hcore
    #             hamming_acc[ind] = hacc
    #             jaccard_core[ind] = jcore
    #             jaccard_acc[ind] = jacc
    #
    #         hamming_core_sims[i] = hamming_core
    #         hamming_acc_sims[i] = hamming_acc
    #         jaccard_core_sims[i] = jaccard_core
    #         jaccard_acc_sims[i] = jaccard_acc

    # determine site specific mutation rates. Get from gamma distribution or manual entries
    gamma_shape = 20.0
    gamma_scale = 1.0

    for i in range(num_sim):

        print("Simulation " + str(i + 1))
        core_ref = np.random.choice(base_choices, size_core, p=base_freq)
        acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)

        # pull out variable sites in core_ref
        sites = np.array(random.sample(range(core_ref.size), k=core_num_var))
        present = np.full(core_ref.size, False)
        present[sites] = True
        core_var = core_ref[present]
        core_invar = core_ref[np.invert(present)]

        hamming_core = [None] * len(core_mu)
        hamming_acc = [None] * len(acc_mu)
        jaccard_core = [None] * len(core_mu)
        jaccard_acc = [None] * len(acc_mu)

        hamming_core = [None] * len(core_mu)
        hamming_acc = [None] * len(acc_mu)
        jaccard_core = [None] * len(core_mu)
        jaccard_acc = [None] * len(acc_mu)

        if core_sites_man is None:
            core_site_mu = calc_gamma(core_num_var, core_sites, gamma_shape, gamma_scale)
        else:
            core_site_mu = calc_man(core_num_var, core_sites_man)

        if acc_sites_man is None:
            acc_site_mu = calc_gamma(size_acc, acc_sites, gamma_shape, gamma_scale)
        else:
            acc_site_mu = calc_man(size_acc, acc_sites_man)

        for ind, val in enumerate(core_mu):
            ind, hcore, hacc, jcore, jacc, acc_vs_core, pangenome_frac = gen_distances(ind, core_var, acc_ref,
                                                                                       core_invar, num_core, core_mu, acc_mu,
                                                                                       adj, avg_gene_freq, base_mu, core_site_mu,
                                                                                       acc_site_mu)
            hamming_core[ind] = hcore
            hamming_acc[ind] = hacc
            jaccard_core[ind] = jcore
            jaccard_acc[ind] = jacc

        hamming_core_sims[i] = hamming_core
        hamming_acc_sims[i] = hamming_acc
        jaccard_core_sims[i] = jaccard_core
        jaccard_acc_sims[i] = jaccard_acc
#
#     mu_rates = (core_mu, acc_mu, core_mu, acc_mu)
#     distances = (hamming_core_sims, hamming_acc_sims, jaccard_core_sims, jaccard_acc_sims)
#     mu_names = ("core_mu", "acc_mu", "core_mu", "acc_mu")
#     distance_names = ("hamming_core", "hamming_acc", "jaccard_core", "jaccard_acc")
#
#     generate_graph(mu_rates, distances, mu_names, distance_names, "./")
#
#     # for var1, var2, name1, name2 in zip((core_mu, acc_mu, core_mu, acc_mu),
#     #                                     (hamming_core, hamming_acc, jaccard_core, jaccard_acc),
#     #                                     ("core_mu", "acc_mu", "core_mu", "acc_mu"),
#     #                                     ("hamming_core", "hamming_acc", "jaccard_core", "jaccard_acc")):
#     #
#     #     #plt.style.use('_mpl-gallery')
#     #
#     #     # make data
#     #     x = np.array(var1)
#     #     y = np.array(var2)
#     #
#     #     # plot
#     #     fig, ax = plt.subplots()
#     #
#     #     ax.plot(x, y, linewidth=2.0)
#     #
#     #     lims = [
#     #         np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#     #         np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
#     #     ]
#     #
#     #     ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
#     #     ax.set_aspect('equal')
#     #     ax.set_xlim(lims)
#     #     ax.set_ylim(lims)
#     #     ax.set_xlabel(name1)
#     #     ax.set_ylabel(name2)
#     #
#     #     fig.savefig(name2 + ".png")
#     #
#     #     plt.clf