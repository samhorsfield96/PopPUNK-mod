
import numpy as np
import argparse
from scipy.spatial.distance import jensenshannon, canberra, cosine, euclidean
from scipy.special import rel_entr, kl_div
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

def get_options():
    description = 'Generate KS distribution for an input file. If two files provided, calculates distance.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python KDE_distance.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile1',
                    required=True,
                    help='Input distance file 1.')
    IO.add_argument('--infile2',
                    default=None,
                    help='Input distance file 2.')
    IO.add_argument('--outpref',
                    default=None,
                    help='Output prefix.')
    
    return parser.parse_args()

def generate_KS_dist(x):
    x_sorted = np.sort(x)
    cdf = np.arange(1, len(x_sorted) + 1) / len(x_sorted)
    return x_sorted, cdf

def KS_dist_2d(df):
    core_sorted, core_cdf = generate_KS_dist(data[:, 0])
    acc_sorted, acc_cdf = generate_KS_dist(data[:, 1])

    return core_sorted, core_cdf, acc_sorted, acc_cdf
    import matplotlib.pyplot as plt

def compute_ecdf(data):
    """Compute empirical CDF of a 1D array."""
    x = np.sort(data)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

def plot_feature_ecdfs(data, label=None, outpref=None):
    n_features = data.shape[1]
    for i in range(n_features):
        x, y = compute_ecdf(data[:, i])
        plt.plot(x, y, label=f'{label} Feature {i}' if label else f'Feature {i}')
    plt.xlabel('Value')
    plt.ylabel('CDF')
    if label:
        plt.title(f"CDFs for {label}")
    plt.legend()
    plt.grid(True)
    if outpref:
        plt.savefig(outpref + ".png", bbox_inches='tight', dpi=300)
    plt.close()

def ks_distances(A, B):
    ks_results = []
    for i in range(A.shape[1]):
        stat, pval = ks_2samp(A[:, i], B[:, i])
        ks_results.append((i, stat, pval))
    return ks_results


def main():
    options = get_options()
    infile1 = options.infile1
    infile2 = options.infile2
    outpref = options.outpref

    df1 = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    df2 = np.loadtxt(infile2, delimiter='\t', dtype='float64')

    if outpref != None:
        outpref1 = outpref + "_df1"
        outpref2 = outpref + "_df2"
    else:
        outpref1, outpref2 = None, None

    plot_feature_ecdfs(df1, label="df1", outpref=outpref1)
    plot_feature_ecdfs(df2, label="df2", outpref=outpref2)

    results = ks_distances(df1, df2)

    # Print results
    for i, stat, pval in results:
        print(f"Feature {i}: KS Statistic = {stat:.3f}, p-value = {pval:.3g}")
   

if __name__ == "__main__":
    main()