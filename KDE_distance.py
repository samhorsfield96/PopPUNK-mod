
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

def js_distance_integral_kde(kde_p, kde_q, eps=1e-12):
    """
    Compute Jensen-Shannon distance between two sklearn KDEs over [0,1]^2 via integral form.
    
    Parameters:
        kde_p, kde_q: fitted sklearn.neighbors.KernelDensity instances
        eps: small value to avoid log(0)
    
    Returns:
        js_distance: float in [0, 1]
    """

    def eval_kde(kde, x, y):
        sample = np.array([[x, y]])
        log_density = kde.score_samples(sample)[0]
        return max(np.exp(log_density), eps)

    def integrand(y, x):
        p = eval_kde(kde_p, x, y)
        q = eval_kde(kde_q, x, y)
        m = max(0.5 * (p + q), eps)
        return 0.5 * (
            p * np.log(2 * p / m) +
            q * np.log(2 * q / m)
        )

    bounds = [(0, 1), (0, 1)]
    js_div, _ = nquad(integrand, bounds)
    js_dist = np.sqrt(js_div)
    return js_dist


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
    
    return parser.parse_args()

def get_kde(X):
    # KDE estimate
    kde = KernelDensity(bandwidth=0.03, metric='euclidean',
                        kernel='epanechnikov', algorithm='ball_tree')
    # convert nan values to 0.0
    X[np.isnan(X)] = 0.0
    #scale = np.amax(X, axis = 0)

    #X /= scale

    # convert nan values to 0.0 again incase of 0 division error
    #X[np.isnan(X)] = 0.0
    
    kde.fit(X)
    return kde

def generate_samples(grid_params, kde):
    xx, yy, xy = grid_params
    
    #print(f"kde.score_samples(xy): {kde.score_samples(xy)}")
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
    # print(f"np.vstack([df1, df2]) {np.vstack([df1, df2])}")

    # print(f"np.max(df1): {np.max(df1, axis=0)}")
    # print(f"np.min(df1): {np.min(df1, axis=0)}")
    # print(f"np.max(df2): {np.max(df2, axis=0)}")
    # print(f"np.min(df2): {np.min(df2, axis=0)}")

    # scale dataframes between 0-1
    df1 = scaler.transform(df1)
    df2 = scaler.transform(df2)

    # print(f"df1: {df1}")
    # print(f"df2: {df2}")

    # print(f"np.max(df1): {np.max(df1, axis=0)}")
    # print(f"np.min(df1): {np.min(df1, axis=0)}")
    # print(f"np.max(df2): {np.max(df2, axis=0)}")
    # print(f"np.min(df2): {np.min(df2, axis=0)}")

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

    df1 = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    df2 = np.loadtxt(infile2, delimiter='\t', dtype='float64')

    # js_distance = KDE_JS_divergence(df1, df2, gamma=0.5, eps=1e-12, log=True)
    # print(f"js_distance_logged: {js_distance}")
    # print(f"log_js_distance_logged: {np.log(js_distance)}")
    #euc_dist = euclidean(np.array([0.0]), np.array([js_distance]))
    #print(f"js_euc_distance: {js_distance}")

    # js_distance = KDE_JS_divergence(df1, df2, gamma=0.5, eps=1e-12, log=False)
    # print(f"js_distance: {js_distance}")
    # print(f"js_distance: {np.log(js_distance)}")

    js_distance = KDE_JS_divergence(df1, df2, gamma=0.25, eps=0.0, log=False, outplot=infile2)
    print(f"js_distance_no_eps: {js_distance}")
    #print(f"log_js_distance_no_eps: {np.log(js_distance)}")

    #KL_distance = np.sum(rel_entr(z1, z2))
    # KL_distance = KDE_KL_divergence(df1, df2)
    # print(f"KL_distance: {KL_distance}")
    # print(f"log_KL_distance: {np.log(KL_distance)}")
    #euc_dist = euclidean(np.array([0.0]), np.array([KL_distance]))
    #print(f"KL_euc_dist: {euc_dist}")

    # baseline_df = np.zeros((df1.shape))
    # print(f"baseline_df: {baseline_df}")
    # scaler = MinMaxScaler()
    # scaler.fit(np.vstack([df1, baseline_df]))

    # # scale dataframes between 0-1
    # df1 = scaler.transform(df1)
    # baseline_df = scaler.transform(baseline_df)

    # print(f"df1: {df1}")
    # print(f"baseline_df: {baseline_df}")

    # # generate grid
    # grid_params = get_grid(0, 1, 100)

    # # get KDE samples
    # z1 = create_KDE_dist(df1, grid_params, 0)
    # z2 = create_KDE_dist(df2, grid_params, 0)

    # print(f"z1: {z1}")
    # print(f"z2: {z2}")

    # KL_distance = KDE_KL_divergence(df1, baseline_df)
    # print(f"baseline KL_distance: {KL_distance}")

if __name__ == "__main__":
    main()