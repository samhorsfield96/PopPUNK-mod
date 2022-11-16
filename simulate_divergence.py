import math
import glob
import os
import pandas as pd
from scipy.spatial import distance
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import random
from math import e
from scipy.stats import gamma
from numba import jit

def read_files(in_dir, prefix=""):
    all_files = glob.glob(os.path.join(in_dir, prefix + "*.txt"))

    li = []
    for filename in all_files:
        #file_pref = os.path.splitext(os.path.basename(filename))[0]
        #split_file_pref = file_pref.split("_distances_sample")
        #name, sample = (split_file_pref[0], split_file_pref[1])

        df = pd.read_csv(filename, index_col=None, header=None, sep="\t")

        # rename columns
        df.rename(columns={df.columns[0]: "Sample1", df.columns[1] : "Sample2", df.columns[2]: "Core",
                           df.columns[3]: "Accessory"}, inplace=True)

        # drop first two columns
        #df['Species'] = name
        #df['Sample'] = int(sample)
        df['Core'] = pd.to_numeric(df['Core'])
        df['Accessory'] = pd.to_numeric(df['Accessory'])

        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)

    return frame

def recurse_prob(x, weight):
    if x == 1:
        return weight
    else:
        return (1 - recurse_prob(x - 1, weight)) * weight

def calc_man_vec(array_size, vec_size, bin_probs, batch_size):
    no_split = np.shape(bin_probs)[1]

    vec_size = np.repeat([vec_size], batch_size)

    # get bins for pdf
    integer = vec_size // no_split
    mod = np.reshape((vec_size % no_split).astype(np.int64), (-1, 1))

    bins = np.transpose(np.array([integer] * np.shape(bin_probs)[1]).astype(int))

    # add modulus to bins to assign all sites (first is vectorised)
    bins[:, 0] += mod[:, 0]
    #bins[0] += mod[0]

    # assign probabilities to sites
    site_mu = np.zeros((np.shape(bin_probs)[0], array_size))

    scaled_bin_probs = bin_probs / bins

    start = np.zeros(np.shape(bin_probs)[0], dtype=np.int64)
    end = np.copy(bins[:, 0])
    #end = np.copy(bins[0])

    r = np.arange(site_mu.shape[1])

    # for testing
    #non_zero = np.zeros(np.shape(bin_probs)[0])

    for index in range(no_split):
        mask = (start[:, None] <= r) & (end[:, None] > r)
        #mask = (start <= r) & (end > r)

        scaled_bins_index = scaled_bin_probs[:, index]

        for entry in range(scaled_bins_index.size):
            site_mu[entry][mask[entry]] = scaled_bins_index[entry]
        #ite_mu[0][mask] = scaled_bins_index

        if index < np.shape(bin_probs)[1] - 1:
            start += bins[:, index]
            end += bins[:, index + 1]
            # start += bins[index]
            # end += bins[index + 1]

        # for testing
        #non_zero += np.count_nonzero(mask, axis=1)

    #for testing
    # for row in range(site_mu.shape[0]):
    #     sum_sites_mu = np.sum(site_mu[row])

    return site_mu

@jit(nopython=True)
def sim_divergence_vec(ref, mu, core, freq, site_mu, negative_value):
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
        if negative_value[i]:
            continue
        for j in range(mu.shape[1]):
            query[i][j] = ref[i].copy()

            num_sites = num_sites_vec[i][j]
            if num_sites > 0:
                total_sites = 0
                while total_sites < num_sites:
                    to_sample = num_sites - total_sites

                    # pick all sites to be mutated (first is vectorised)
                    sites = index_array[np.searchsorted(np.cumsum(site_mu[i]), np.random.rand(to_sample))]
                    #sites = index_array[np.searchsorted(np.cumsum(site_mu), np.random.rand(to_sample))]

                    unique = np.array(list(set(sites)))
                    counts = np.array([np.count_nonzero(sites == val) for val in unique])

                    max_count = np.max(counts)

                    for count in range(1, max_count + 1):
                        # determine number of times each site can be mutated
                        sample_sites = unique[counts >= count]

                        # determine sites with and without change (first is vectorised)
                        changes = choices[np.searchsorted(np.cumsum(freq[i]), np.random.rand(sample_sites.size))]
                        #changes = choices[np.searchsorted(np.cumsum(freq), np.random.rand(sample_sites.size))]

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

def sim_divergence(ref, mu, core, freq, site_mu, dispersion):
    # add dispersion to number of sites based on normal distribution
    mu = abs(np.random.normal(mu, dispersion))

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

            unique, counts = np.unique(sites, return_counts=True)

            max_count = np.max(counts)

            for count in range(1, max_count + 1):
                # determine number of times each site can be mutated
                sample_sites = unique[counts >= count]

                # determine sites with and without change
                changes = np.random.choice(choices, sample_sites.size, p=freq)

                non_mutated = query[sample_sites] == changes

                query[sample_sites] = changes

                # determine number of actual changes
                total_sites += changes.size - np.count_nonzero(non_mutated)

    return query, total_sites

def gen_distances(index, core_var, acc_var, core_invar, num_core, core_mu, acc_mu, adj, avg_gene_freq, base_mu,
                  core_site_mu, acc_site_mu, sim_core_dispersion, sim_acc_dispersion):
    sites_mutated = []

    # mutate genomes
    core_query1, total_sites = sim_divergence(core_var, core_mu[index], True, base_mu, core_site_mu, sim_core_dispersion)
    sites_mutated.append(total_sites)
    #core_query2, total_sites = sim_divergence(core_var, core_mu[index], True, base_mu, core_site_mu, sim_core_dispersion)
    #sites_mutated.append(total_sites)
    acc_query1, total_sites = sim_divergence(acc_var, acc_mu[index], False, avg_gene_freq, acc_site_mu, sim_acc_dispersion)
    sites_mutated.append(total_sites)
    #acc_query2, total_sites = sim_divergence(acc_var, acc_mu[index], False, avg_gene_freq, acc_site_mu, sim_acc_dispersion)
    #sites_mutated.append(total_sites)

    if adj == True:
        # add core genes to accessory distances
        acc_var = np.append(acc_var, np.ones(num_core))
        acc_query1 = np.append(acc_query1, np.ones(num_core))
        #acc_query2 = np.append(acc_query2, np.ones(num_core))

        # add core invariant sites to core alignments
        core_var = np.append(core_var, core_invar)
        core_query1 = np.append(core_query1, core_invar)
        #core_query2 = np.append(core_query2, core_invar)

    acc_vs_core = 1

    if core_mu[index] != 0:
        # determine accessory vs core divergence rate
        prop_subs_core = (sites_mutated[0]) / (core_query1.size + core_var.size)
        prop_subs_acc = (sites_mutated[1]) / (acc_query1.size + acc_var.size)
        acc_vs_core = prop_subs_acc / prop_subs_core

    # determine pangenome_frac
    match = acc_query1 == acc_var
    zeros_match = match[acc_query1 == 0]
    num_zero_match = np.count_nonzero(zeros_match)

    pangenome_frac = (acc_query1.size - num_zero_match) / acc_query1.size

    hamming_core = distance.hamming(core_query1, core_var)
    hamming_acc = distance.hamming(acc_query1, acc_var)
    jaccard_core = distance.jaccard(core_query1, core_var)
    jaccard_acc = distance.jaccard(acc_query1, acc_var)

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
        c, cov = curve_fit(model, reg_x, reg_y, maxfev=5000,
                           bounds=(([0, 0, 0, -np.inf]), (np.inf, 1, np.inf, 0)))
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

        fig.savefig(outpref + "_" + name + "_vs_pangenome_frac.png")

        ax.clear()
        plt.clf
        plt.cla

def generate_graph(mu_rates, distances, mu_names, distance_names, outpref, core_adj, acc_adj, adjusted):
    for var1, var2, name1, name2 in zip(mu_rates, distances, mu_names, distance_names):
        # plot
        fig, ax = plt.subplots()

        # make data, as comparing two diverged sequences, multiply by 2
        x = np.array(var1)
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
        #ax.set_xlim(lims)
        #ax.set_ylim(0, 1)
        ax.set_xlabel(name1)
        ax.set_ylabel(name2)
        #ax.legend()

        fig.savefig(outpref + "_" + name2 + ".png")

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

        #ax.plot(x, y, linewidth=2.0, label="Sim" + str(sim))
        ax.scatter(x, y)

        reg_x[(sim - 1) * x.size : ((sim - 1) * x.size) + x.size] = x
        reg_y[(sim - 1) * y.size : ((sim - 1) * y.size) + y.size] = y

        sim += 1

    #fit model, determine uncertainty
    # c, cov = curve_fit(model, reg_x, reg_y, maxfev=5000, bounds=(([0,0,0,-np.inf]), (np.inf,1,np.inf,0)))
    # d_c0 = np.sqrt(cov[0][0])
    # d_c1 = np.sqrt(cov[1][1])
    # d_c2 = np.sqrt(cov[2][2])
    # d_c3 = np.sqrt(cov[3][3])

    # predict using new model
    # x = np.array(distances[0][0])
    # y = np.array([model(j, c[0], c[1], c[2], c[3]) for j in x])
    # ax.plot(x, y, linewidth=2.0, label="Model")
    # print("Accessory vs. core model parameters:")
    # print(c)
    # print("[" + str(d_c0) + " " + str(d_c1) + " " + str(d_c2) + " " + str(d_c3) + "]")

    # with open(outpref + "model_parameters.txt", "w") as f:
    #     f.write(np.array2string(c))
    #     f.write("[" + str(d_c0) + " " + str(d_c1) + " " + str(d_c2) + " " + str(d_c3) + "]")

    lims = [
        np.min([ax.get_xlim(), ax.get_xlim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_xlim()]),  # max of both axes
    ]

    # ax.set_ylim(0, 1)
    # ax.set_xlim(lims)
    ax.set_xlabel("Core Hamming distance")
    ax.set_ylabel("Accessory Jaccard Distance")
    #ax.legend()

    fig.savefig(outpref + "_core_vs_acc.png")
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
    #ax.legend()

    fig.savefig(outpref + "_acc_hamming_vs_jaccard.png")
    plt.cla

    plt.clf