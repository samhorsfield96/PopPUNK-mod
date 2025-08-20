
import numpy as np
import argparse
from scipy.spatial.distance import jensenshannon, canberra, cosine, euclidean
from scipy.special import rel_entr, kl_div
from sklearn.preprocessing import MinMaxScaler
from scipy.integrate import nquad
import matplotlib.pyplot as plt

try:  # sklearn >= 0.22
    from sklearn.neighbors import KernelDensity
except ImportError:
    from sklearn.neighbors.kde import KernelDensity


# Copyright John Lees and Nicholas Croucher 2025
def get_grid(minimum, maximum, resolution):
    x = np.linspace(minimum, maximum, resolution)
    y = np.linspace(minimum, maximum, resolution)
    xx, yy = np.meshgrid(x, y)
    xy = np.vstack([yy.ravel(), xx.ravel()]).T

    return(xx, yy, xy)

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
    IO.add_argument('--gamma',
                    default=0.25,
                    type=float,
                    help='Power to raise KDE distributions to for smoothing (Default = 0.25).')
    IO.add_argument('--plot-outpref',
                    default=None,
                    type=str,
                    help='Prefix for output plots. If not provided, no plots generated.')
    
    return parser.parse_args()

def get_kde(X):
    # KDE estimate
    kde = KernelDensity(bandwidth=0.03, metric='euclidean',
                        kernel='epanechnikov', algorithm='ball_tree')
    # convert nan values to 0.0
    X[np.isnan(X)] = 0.0
    
    kde.fit(X)
    return kde

def generate_samples(grid_params, kde):
    xx, yy, xy = grid_params
    
    z = np.exp(kde.score_samples(xy))
    z = z.reshape(xx.shape)

    return z

def create_KDE_dist(df, grid_params, gamma, eps):
    kde = get_kde(df)
    z = generate_samples(grid_params, kde)

    # reweight distribution
    z = np.power(z, gamma)

    # normalise
    z = z.ravel()
    z += eps
    z /= z.sum()

    return z

def scale_KDE(df1, df2, gamma, eps):
    scaler = MinMaxScaler()
    scaler.fit(np.vstack([df1, df2]))

    # scale dataframes between 0-1
    df1 = scaler.transform(df1)
    df2 = scaler.transform(df2)

    # generate grid
    grid_params = get_grid(0, 1, 100)

    # get KDE samples
    z1 = create_KDE_dist(df1, grid_params, gamma, eps)
    z2 = create_KDE_dist(df2, grid_params, gamma, eps)

    return z1, z2, grid_params

def KDE_KL_divergence(df1, df2, gamma=1.0, eps=1e-12):
    z1, z2 = scale_KDE(df1, df2, gamma, eps)

    # calculate KL divergence
    KL_divergence = np.sum(rel_entr(z1, z2))
    return KL_divergence

def KDE_JS_divergence(df1, df2, gamma=1.0, eps=1e-12, log=False, outplot=None):
    if log:
        df1, df2 = np.log(df1 + eps), np.log(df2 + eps)
    z1, z2, grid_params = scale_KDE(df1, df2, gamma, eps)

    # plot contours
    if outplot != None:
        scatter_alpha = 1
        z1_square = z1.reshape(grid_params[0].shape).T
        plt.figure(figsize=(11, 8), dpi= 160, facecolor='w', edgecolor='k')
        levels = np.linspace(z1_square.min(), z1_square.max(), 100)
        plt.contour(grid_params[0], grid_params[1], z1_square, levels=levels[1:], cmap='plasma')
        plt.scatter(df1[:,0].flat, df1[:,1].flat, s=1, alpha=scatter_alpha)
        plt.savefig(outplot + '_z1_contours.png')
        plt.close()

        z2_square = z2.reshape(grid_params[0].shape).T
        plt.figure(figsize=(11, 8), dpi= 160, facecolor='w', edgecolor='k')
        levels = np.linspace(z2_square.min(), z2_square.max(), 100)
        plt.contour(grid_params[0], grid_params[1], z2_square, levels=levels[1:], cmap='plasma')
        plt.scatter(df2[:,0].flat, df2[:,1].flat, s=1, alpha=scatter_alpha)
        plt.savefig(outplot + '_z2_contours.png')
        plt.close()

    # calculate KL divergence
    js_distance = jensenshannon(z1, z2, axis=0)
    return js_distance

def main():
    options = get_options()
    infile1 = options.infile1
    infile2 = options.infile2
    gamma = options.gamma
    plot_outpref = options.plot_outpref

    df1 = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    df2 = np.loadtxt(infile2, delimiter='\t', dtype='float64')

    
    js_distance = KDE_JS_divergence(df1, df2, gamma=gamma, eps=0.0, log=False, outplot=plot_outpref)
    print(f"js_distance_no_eps: {js_distance}")


if __name__ == "__main__":
    main()