import argparse

import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.cluster import HDBSCAN
from sklearn import utils
from collections import defaultdict
from pathlib import Path
import seaborn as sns
import random
from multiprocessing import Pool
from functools import partial

try:  # sklearn >= 0.22
    from sklearn.neighbors import KernelDensity
except ImportError:
    from sklearn.neighbors.kde import KernelDensity

def get_options():
    description = 'Compare parameters from a gridsearch of Pansim runs'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python compare_gridsearch.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    required=True,
                    help='Directory containing run files from gridsearch run')
    IO.add_argument('--outpref',
                    default = "output",
                    help='Output prefix. Default = "output"')
    IO.add_argument('--plot',
                    action="store_true",
                    default=False,
                    help='Generate plots.')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1.')

    return parser.parse_args()

# Copyright John Lees and Nicholas Croucher 2025
def get_grid(minimum, maximum, resolution):
    x = np.linspace(minimum, maximum, resolution)
    y = np.linspace(minimum, maximum, resolution)
    xx, yy = np.meshgrid(x, y)
    xy = np.vstack([yy.ravel(), xx.ravel()]).T

    return(xx, yy, xy)

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
    plt.savefig(out_prefix + '_distanceDistribution.png')
    plt.close()

# fit asymptotic curve using exponential decay
def negative_exponential(x, b0, b1, b2): # based on https://isem-cueb-ztian.github.io/Intro-Econometrics-2017/handouts/lecture_notes/lecture_10/lecture_10.pdf and https://www.statforbiology.com/articles/usefulequations/
    return b0 * (1 - np.exp(-b1 * x)) + b2

def generate_summary_stats(file, plot, outpref):
    filename = os.path.splitext(os.path.basename(file))[0]

    split_filename = filename.split("_")

    # get all variables for file:
    values_dict = {}
    values_dict["competition"] = False
    for i in range(len(split_filename) - 2):
        current_index = split_filename[i]

        try:
            next_index = float(split_filename[i + 1])
            values_dict[current_index] = next_index
        except:
            continue

    if "non_competition" not in filename and "competition" in filename:
        values_dict["competition"] = True
    
    obs_df = np.loadtxt(file, delimiter='\t', dtype='float64')

    values_dict["b0"] = "NA"
    values_dict["b1"] = "NA"
    values_dict["b2"] = "NA"
    values_dict["b0_err"] = "NA"
    values_dict["b1_err"] = "NA"
    values_dict["b2_err"] = "NA"

    values_dict["max_core"] = np.max(obs_df[:,0])
    values_dict["min_core"] = np.min(obs_df[:,0])
    values_dict["max_acc"] = np.max(obs_df[:,1])
    values_dict["min_acc"] = np.min(obs_df[:,1])
    values_dict["mean_acc"] = np.mean(obs_df[:,1])
    values_dict["mean_core"] = np.mean(obs_df[:,0])
    
    values_dict["median_acc"] = np.median(obs_df[:,1])
    values_dict["median_core"] = np.median(obs_df[:,0])
    values_dict["quant_25_acc"] = np.quantile(obs_df[:,1], 0.25)
    values_dict["quant_25_core"] = np.quantile(obs_df[:,0], 0.25)
    values_dict["quant_75_acc"] = np.quantile(obs_df[:,1], 0.75)
    values_dict["quant_75_core"] = np.quantile(obs_df[:,0], 0.75)

    # cluster positions
    values_dict["min_centroid_core"] = "NA"
    values_dict["min_centroid_acc"] = "NA"
    values_dict["max_centroid_acc"] = "NA"
    values_dict["max_centroid_core"] = "NA"
    values_dict["quant_25_centroid_core"] = "NA"
    values_dict["quant_25_centroid_acc"] = "NA"
    values_dict["quant_75_centroid_core"] = "NA"
    values_dict["quant_75_centroid_acc"] = "NA"
    values_dict["median_centroid_core"] = "NA"
    values_dict["median_centroid_acc"] = "NA"

    # cluster sizes
    values_dict["min_centroid_size"] = "NA"
    values_dict["max_centroid_size"] = "NA"
    values_dict["quant_25_centroid_size"] = "NA"
    values_dict["quant_75_centroid_size"] = "NA"
    values_dict["median_centroid_size"] = "NA"

    # try negative_exponential model fit
    try:
        popt, pcov = curve_fit(negative_exponential, obs_df[:,0], obs_df[:,1], p0=[1.0, 1.0, 0.0], bounds=([0.0, 0.0, 0.0], [1.0, np.inf, 1.0]))
        b0, b1, b2 = popt
        # get 1 std deviation error of parameters
        b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))

        values_dict["b0"] = b0
        values_dict["b1"] = b1
        values_dict["b2"] = b2
        values_dict["b0_err"] = b0_err
        values_dict["b1_err"] = b1_err
        values_dict["b2_err"] = b2_err

        x_fit = np.linspace(0, obs_df[:,0].max(), 100)
        y_fit = negative_exponential(x_fit, *popt)

        if plot:
            plot_scatter(obs_df, outpref + "_" + filename, x_fit, y_fit)
    except:
        pass
    
    # try HDSCAN, absolute minimum 10 samples per cluster
    min_cluster_size = max(round(len(obs_df[:,1]) * 0.01), 10)
    min_samples = max(round(len(obs_df[:,1]) * 0.0001), 10)
    try: 
        hdb = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, algorithm="ball_tree")
        hdb.fit(obs_df)

        # Get unique cluster labels (excluding noise, which is labeled as -1)
        labels = hdb.labels_
        unique_labels = np.unique(labels[labels != -1])  # Ignore noise (-1)

        # Compute cluster sizes and centroids
        data = []

        for label in unique_labels:
            hdb_points = obs_df[labels == label]
            centroid = hdb_points.mean(axis=0)
            size = len(hdb_points)
            data.append([size, centroid[0], centroid[1]])

        # Convert to NumPy array
        hdb_data = np.array(data)

        # Sort by cluster size
        hdb_data = hdb_data[hdb_data[:, 0].argsort()]  # Sort by first column (size)
        

        # Extract sizes and centroids
        sizes = hdb_data[:, 0]       # Cluster sizes
        centroids = hdb_data[:, 1:]  # Centroid positions (x, y)

        # Extract min/max cluster size and corresponding centroids
        min_size, min_centroid = sizes[0], centroids[0]
        max_size, max_centroid = sizes[-1], centroids[-1]

        # Compute percentiles
        size_percentiles = np.percentile(sizes, [25, 50, 75])
        centroid_percentiles = np.percentile(centroids, [25, 50, 75], axis=0)

        # Extract percentile values
        size_25th, size_median, size_75th = size_percentiles
        centroid_25th, centroid_median, centroid_75th = centroid_percentiles

        # cluster positions
        values_dict["min_centroid_core"] = min_centroid[0]
        values_dict["min_centroid_acc"] = min_centroid[1]
        values_dict["max_centroid_acc"] = max_centroid[0]
        values_dict["max_centroid_core"] = max_centroid[1]
        values_dict["quant_25_centroid_core"] = centroid_25th[0]
        values_dict["quant_25_centroid_acc"] = centroid_25th[1]
        values_dict["quant_75_centroid_core"] = centroid_75th[0]
        values_dict["quant_75_centroid_acc"] = centroid_75th[1]
        values_dict["median_centroid_core"] = centroid_median[0] 
        values_dict["median_centroid_acc"] = centroid_median[1]

        # cluster sizes
        values_dict["min_centroid_size"] = int(min_size) / len(obs_df[:,1])
        values_dict["max_centroid_size"] = int(max_size) / len(obs_df[:,1])
        values_dict["quant_25_centroid_size"] = int(size_25th) / len(obs_df[:,1])
        values_dict["quant_75_centroid_size"] = int(size_75th) / len(obs_df[:,1])
        values_dict["median_centroid_size"] = int(size_median) / len(obs_df[:,1])

        #print(values_dict)

        # Plot clustered points
        if plot:
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x=obs_df[:, 0], y=obs_df[:, 1], hue=labels, palette="viridis", s=5, alpha=0.6, edgecolor=None)

            # Mark centroids
            centroid_markers = {
                "Smallest": ("grey", "o"),
                "Largest": ("grey", "p"),
                "25th %": ("grey", "s"),
                "Median": ("grey", "D"),
                "75th %": ("grey", "^"),
            }

            for name, (color, marker), centroid in zip(
                centroid_markers.keys(),
                centroid_markers.values(),
                [min_centroid, max_centroid, centroid_25th, centroid_median, centroid_75th]
            ):
                plt.scatter(*centroid, color=color, marker=marker, s=200, label=name, edgecolor="black")

            # Customize plot
            plt.legend(title="Centroids", loc="best")
            plt.title("HDBSCAN Clustering with Centroid Markers")
            plt.xlabel("Feature 1")
            plt.ylabel("Feature 2")
            plt.grid(True)

            # Show plot
            plt.savefig(outpref + "_" + filename + "_HBSCAN.png")
            plt.close()
    except:
        pass
            
    return values_dict

def main():
    options = get_options()
    indir = options.indir
    outpref = options.outpref
    plot = options.plot
    threads = options.threads

    files = []
    for path in Path(indir).glob("**/*.tsv"):
        # Print the path (file or directory) to the console
        files.append(str(path))
    
    # print(indir)
    # print(files)

    #print(files)
    
    value_dict_list = []
    with Pool(processes=threads) as pool:
        for result_dict in pool.map(partial(generate_summary_stats, plot=plot, outpref=outpref), files):
            value_dict_list.append(result_dict)
    
    df = pd.DataFrame.from_dict(value_dict_list)
    #print(df)
    df.to_csv(outpref + "_results.tsv", sep = "\t", index = False)


if __name__ == "__main__":
    main()