import pandas as pd
import glob
import os
from multiprocessing import Pool
from functools import partial
from simulate_divergence import *
import tqdm
from sklearn.metrics import mean_squared_error
import math

def read_files(in_dir, prefix=""):
    all_files = glob.glob(os.path.join(in_dir, "*.txt"))

    li = []
    for filename in all_files:
        file_pref = os.path.splitext(os.path.basename(filename))[0]
        split_file_pref = file_pref.split("_distances_sample")
        name, sample = (split_file_pref[0], split_file_pref[1])

        df = pd.read_csv(filename, index_col=None, header=None, sep="\t")

        # rename columns
        df.columns = ['Species', 'Sample', 'Core', 'Accessory']

        # drop first two columns
        df['Species'] = name
        df['Sample'] = sample

        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)

    return frame

if __name__ == "__main__":
    df = read_files("/mnt/c/Users/sth19/PycharmProjects/PhD_project/distance_sim/distances")
    outpref = "./"

    threads = 4
    size_core = 1140000
    num_core = 1194
    num_pangenome = 5442
    num_sim = 2
    core_invar = 106196
    adjusted = True

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

    adj = True
    core_num_var = 106196

    hamming_core = [None] * len(core_mu)
    hamming_acc = [None] * len(core_mu)
    jaccard_core = [None] * len(core_mu)
    jaccard_acc = [None] * len(core_mu)

    # generate references
    core_ref = np.random.choice(base_choices, size_core, p=base_freq)
    size_acc = num_pangenome - num_core
    acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)


    # specify range of accesory vs. core rates
    a_vs_c_rates = [0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6, 8, 10]
    #a_vs_c_rates = [1]

    # hold coefficients
    curve_coefficents = [None] * len(a_vs_c_rates)

    for index, rate in enumerate(a_vs_c_rates):
        print("Accessory vs. core rate: " + str(rate))
        # determine acc_mu as a function of core_mu
        acc_mu = [x * rate for x in core_mu]

        hamming_core_sims = [None] * num_sim
        hamming_acc_sims = [None] * num_sim
        jaccard_core_sims = [None] * num_sim
        jaccard_acc_sims = [None] * num_sim

        for i in range(num_sim):
            print("Simulation " + str(i + 1))
            core_ref = np.random.choice(base_choices, size_core, p=base_freq)
            acc_ref = np.random.choice(gene_choices, size_acc, p=gene_freq)

            # pull out variable sites in core_ref
            sites = np.array(random.choices(range(core_ref.size), k=core_num_var))
            present = np.full(core_ref.size, False)
            present[sites] = True
            core_var = core_ref[present]
            core_invar = core_ref[np.invert(present)]

            hamming_core = [None] * len(core_mu)
            hamming_acc = [None] * len(acc_mu)
            jaccard_core = [None] * len(core_mu)
            jaccard_acc = [None] * len(acc_mu)

            with Pool(processes=threads) as pool:
                for ind, hcore, hacc, jcore, jacc in tqdm.tqdm(pool.imap(
                        partial(gen_distances, core_var=core_var, acc_ref=acc_ref, core_invar=core_invar,
                                num_core=num_core, core_mu=core_mu, acc_mu=acc_mu, adj=adjusted, avg_gene_freq=avg_gene_freq,
                                base_mu=base_mu),
                        range(0, len(core_mu))), total=len(core_mu)):
                    hamming_core[ind] = hcore
                    hamming_acc[ind] = hacc
                    jaccard_core[ind] = jcore
                    jaccard_acc[ind] = jacc

            hamming_core_sims[i] = hamming_core
            hamming_acc_sims[i] = hamming_acc
            jaccard_core_sims[i] = jaccard_core
            jaccard_acc_sims[i] = jaccard_acc

        # generate a fit for model
        c = fit_cvsa_curve(hamming_core_sims, jaccard_acc_sims)

        curve_coefficents[index] = c

    # create vectorised version of model
    model_vec = np.vectorize(model)

    # create numpy arrays of inputs
    true_x = df['Core'].to_numpy()
    true_y = df['Accessory'].to_numpy()

    # fit to real data, calculate RMSE
    RMSE_results = [None] * len(a_vs_c_rates)
    for index, entry in enumerate(curve_coefficents):
        c = entry
        pred_y = model_vec(true_x, c[0], c[1])
        MSE = mean_squared_error(true_y, pred_y)
        RMSE_results[index] = math.sqrt(MSE)

    # overlap plots
    fig, ax = plt.subplots()
    x = true_x
    y = true_y
    ax.scatter(x, y, alpha=0.15)

    for index in range(len(a_vs_c_rates)):
        c = curve_coefficents[index]
        rate = a_vs_c_rates[index]
        print("Rate: " + str(rate))
        print("Model parameters:\n" + "Accessory vs. core rate " + str(c[0]) + "\n01s, 10s and 11s over pangenome size: " + str(c[1]))
        print("RMSE: " + str(RMSE_results[index]))

        # sample 10 even point across 0 and max core value
        step = np.max(true_x) / 10
        start = 0
        x = []
        for i in range(11):
            x.append(i * step)

        x = np.array(x)
        y = model_vec(x, c[0], c[1])

        ax.plot(x, y, linewidth=2.0, label="Rate: " + str(rate))

    ax.set_xlabel("Core")
    ax.set_ylabel("Accessory")
    ax.legend()

    fig.savefig(outpref + "core_vs_accesory_fit.png")

    plt.clf

    test = 1