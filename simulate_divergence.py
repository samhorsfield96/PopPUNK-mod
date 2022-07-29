from scipy.spatial import distance
import matplotlib.pyplot as plt
import numpy as np
import random
from math import e


def sim_divergence(query, mu, num_letters, core):
    choices = []
    # generate set of choices of alleles to choose from, dropping single allele each time
    if core:
        for k in range(1, num_letters + 1):
            choices.append(np.array([v for v in range(1, num_letters + 1) if v != k]))
    else:
        for k in range(0, num_letters):
            choices.append(np.array([v for v in range(0, num_letters) if v != k]))

    num_sites = round(len(query) * mu)

    #pick all sites to be mutated
    sites = np.array(random.choices(range(query.size), k=num_sites))

    # determine number of times each site can be mutated
    unique, counts = np.unique(sites, return_counts=True)

    # determine greatest number of sites re-mutated
    max_counts = 1
    if sites.size > 0:
        max_counts = np.max(np.unique(counts))

        # iterate over
        for i in range(1, max_counts + 1):
            count_sites = counts >= i
            sample_sites = unique[count_sites]

            # randomly assign base based on element at that position
            # core is 1 indexed, covert to 0 indexed
            if core:
                changes = np.array([np.random.choice(choices[int(query[j]) - 1]) for j in sample_sites])
            else:
                changes = np.array([np.random.choice(choices[int(query[j])]) for j in sample_sites])
            query[sample_sites] = changes

    return query

def gen_distances(index, core_var, acc_ref, core_invar, num_core, core_mu, acc_mu, adj):
    # mutate genomes
    core_query = sim_divergence(np.copy(core_var), core_mu[index], 4, True)
    acc_query = sim_divergence(np.copy(acc_ref), acc_mu[index], 2, False)

    if adj == True:
        # add core genes to accessory distances
        acc_ref = np.append(acc_ref, np.ones(num_core))
        acc_query = np.append(acc_query, np.ones(num_core))

        # add core invariant sites to core alignments
        core_var = np.append(core_var, core_invar)
        core_query = np.append(core_query, core_invar)

    hamming_core = distance.hamming(core_var, core_query)
    hamming_acc = distance.hamming(acc_ref, acc_query)
    jaccard_core = distance.jaccard(core_var, core_query)
    jaccard_acc = distance.jaccard(acc_ref, acc_query)

    return (index, hamming_core, hamming_acc, jaccard_core, jaccard_acc)


def generate_graph(mu_rates, distances, mu_names, distance_names, lengths, outpref, core_adj, adjusted):
    for var1, var2, name1, name2, length in zip(mu_rates, distances, mu_names, distance_names, lengths):

        #plt.style.use('_mpl-gallery')

        # plot
        fig, ax = plt.subplots()

        # make data
        x = np.array(var1)
        for j in range(len(var2)):
            y = np.array(var2[j])

            ax.plot(x, y, linewidth=2.0, label="Sim" + str(j + 1))

        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]

        # plot Jukes-cantor relationship
        if "core" in name1:
            if adjusted:
                y = (3/4 * (1 - e ** (-(4/3) * (x)))) / core_adj
            else:
                y = 3 / 4 * (1 - e ** (-(4 / 3) * (x)))
        else:
            y = 1/2 * (1 - e ** (-(2/1) * x))

        ax.plot(x, y, linewidth=2.0, label="Jukes-Cantor")

        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0, label="y=x")
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(0, 1)
        ax.set_xlabel(name1)
        ax.set_ylabel(name2)
        ax.legend()

        fig.savefig(outpref + name2 + ".png")

        plt.clf

# if __name__ == "__main__":
#     threads = 4
#     size_core = 1140000
#     num_core = 1194
#     num_pangenome = 5442
#     num_sim = 2
#     core_invar = 106196
#
#     # base frequencies are alphabetical, A, C, G, T
#     base_freq = [0.25, 0.25, 0.25, 0.25]
#     base_choices = [1, 2, 3, 4]
#
#     # gene presence/absence frequencies are in order 0, 1
#     gene_freq = [0.5, 0.5]
#     gene_choices = [0, 1]
#
#     # core mu is number of differences per base of alignment
#     core_mu = [0.1 * i for i in range(0, 21, 2)]
#
#     # acc mu is number of differences per gene in accessory genome
#     acc_mu = [0.1 * i for i in range(0, 11)]
#
#     hamming_core = [None] * len(core_mu)
#     hamming_acc = [None] * len(core_mu)
#     jaccard_core = [None] * len(core_mu)
#     jaccard_acc = [None] * len(core_mu)
#
#     # generate references
#     core_ref = np.random.choice(base_choices, size_core, p=base_freq)
#     size_acc = num_pangenome - num_core
#     acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)
#
#     hamming_core_sims = [None] * num_sim
#     hamming_acc_sims = [None] * num_sim
#     jaccard_core_sims = [None] * num_sim
#     jaccard_acc_sims = [None] * num_sim
#
#     with Pool(processes=threads) as pool:
#         for i in range(num_sim):
#             print("Simulation " + str(i + 1))
#
#             hamming_core = [None] * len(core_mu)
#             hamming_acc = [None] * len(acc_mu)
#             jaccard_core = [None] * len(core_mu)
#             jaccard_acc = [None] * len(acc_mu)
#
#
#             for ind, hcore, hacc, jcore, jacc in tqdm.tqdm(pool.imap(
#                     partial(gen_distances, core_ref=core_ref, acc_ref=acc_ref, core_invar=core_invar,
#                             num_core=num_core, core_mu=core_mu, acc_mu=acc_mu),
#                     range(0, len(core_mu))), total=len(core_mu)):
#                 hamming_core[ind] = hcore
#                 hamming_acc[ind] = hacc
#                 jaccard_core[ind] = jcore
#                 jaccard_acc[ind] = jacc
#
#             hamming_core_sims[i] = hamming_core
#             hamming_acc_sims[i] = hamming_acc
#             jaccard_core_sims[i] = jaccard_core
#             jaccard_acc_sims[i] = jaccard_acc
#
#     # for ind, val in enumerate(core_mu):
#     #     ind, hcore, hacc, jcore, jacc = gen_distances(ind, core_ref, acc_ref, core_invar, num_core, core_mu, acc_mu)
#     #     hamming_core[ind] = hcore
#     #     hamming_acc[ind] = hacc
#     #     jaccard_core[ind] = jcore
#     #     jaccard_acc[ind] = jacc
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