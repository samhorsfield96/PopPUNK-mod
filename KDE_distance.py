
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

def create_KDE_dist(df, grid_params):
    kde = get_kde(df)
    z = generate_samples(grid_params, kde)
    z = z.ravel()
    z /= z.sum()

    return z

def main():
    options = get_options()
    infile1 = options.infile1
    infile2 = options.infile2

    infile1_df = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    infile2_df = np.loadtxt(infile2, delimiter='\t', dtype='float64')

    grid_params = get_grid(0, 1, 100)

    z1 = create_KDE_dist(infile1_df, grid_params)
    z2 = create_KDE_dist(infile2_df, grid_params)

    js_distance = jensenshannon(z1, z2, axis=0)
    print(f"js_distance: {js_distance}")

if __name__ == "__main__":
    main()