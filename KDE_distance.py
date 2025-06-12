
import numpy as np
import argparse
from run_elfi import get_grid
from scipy.spatial.distance import jensenshannon
from scipy.special import rel_entr

try:  # sklearn >= 0.22
    from sklearn.neighbors import KernelDensity
except ImportError:
    from sklearn.neighbors.kde import KernelDensity

def get_options():
    description = 'Fit KDE to two 2D distributions and generates a distance.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python KDE_distance.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile1',
                    required=True,
                    help='Input distance file 1.')
    IO.add_argument('--infile2',
                    required=True,
                    help='Input distance file 2.')
    
    return parser.parse_args()

def get_kde(X):
    # KDE estimate
    kde = KernelDensity(bandwidth=0.03, metric='euclidean',
                        kernel='epanechnikov', algorithm='ball_tree')
    kde.fit(X)
    return kde

def generate_samples(grid_params, kde):
    xx, yy, xy = grid_params
    
    z = np.exp(kde.score_samples(xy))
    z = z.reshape(xx.shape).T

    return z


def main():
    options = get_options()
    infile1 = options.infile1
    infile2 = options.infile2

    infile1_df = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    infile2_df = np.loadtxt(infile2, delimiter='\t', dtype='float64')

    kde1 = get_kde(infile1_df)
    kde2 = get_kde(infile2_df)

    grid_params = get_grid(0, 1, 100)
    xx, yy, xy = grid_params

    print(f"xx: {xx}")
    print(f"yy: {yy}")
    print(f"xy: {xy}")

    z1 = generate_samples(grid_params, kde1)
    z2 = generate_samples(grid_params, kde2)

    z1 /= z1.sum()
    z2 /= z2.sum()

    print(f"z1: {z1}")
    print(f"z2: {z2}")

    # Compute KL divergence
    kl_pq = np.sum(rel_entr(z1, z2))
    print(f"kl_pq: {kl_pq}")

    # Compute Jensen-Shannon divergence (symmetric and bounded)
    js_divergence = 0.5 * np.sum(rel_entr(z1, 0.5*(z1+z2))) + 0.5 * np.sum(rel_entr(z2, 0.5*(z1+z2)))
    print(f"js_divergence: {js_divergence}")

    # Convert to Jensen-Shannon distance (metric)
    js_distance = np.sqrt(js_divergence)
    print(f"js_distance: {js_distance}")

if __name__ == "__main__":
    main()