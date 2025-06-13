
import numpy as np
import argparse
from scipy.spatial.distance import jensenshannon, canberra, cosine, euclidean
from scipy.special import rel_entr, kl_div
from sklearn.preprocessing import MinMaxScaler
from scipy.integrate import nquad

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
    kde.fit(X)
    return kde

def generate_samples(grid_params, kde):
    xx, yy, xy = grid_params
    
    #print(f"kde.score_samples(xy): {kde.score_samples(xy)}")
    z = np.exp(kde.score_samples(xy))
    z = z.reshape(xx.shape).T

    return z

def create_KDE_dist(df, grid_params, eps=1e-12):
    scale = np.amax(df, axis = 0)
    df /= scale

    kde = get_kde(df)
    z = generate_samples(grid_params, kde)
    #print(f"z: {z}")
    z = z.ravel()
    z += eps
    z /= z.sum()

    return z

def scale_KDE(df1, df2, eps=1e-12):
    scaler = MinMaxScaler()
    scaler.fit(np.vstack([df1, df2]))

    # scale dataframes between 0-1
    df1 = scaler.transform(df1)
    df2 = scaler.transform(df2)

    # generate grid
    grid_params = get_grid(0, 1, 100)

    # get KDE samples
    z1 = create_KDE_dist(df1, grid_params, eps)
    z2 = create_KDE_dist(df2, grid_params, eps)

    return z1, z2

def KDE_KL_divergence(df1, df2, eps=1e-12):
    z1, z2 = scale_KDE(df1, df2, eps)

    # calculate KL divergence
    KL_divergence = np.sum(rel_entr(z1, z2))
    return KL_divergence

def KDE_JS_divergence(df1, df2, eps=1e-12):
    z1, z2 = scale_KDE(df1, df2, eps)

    # calculate KL divergence
    js_distance = jensenshannon(z1, z2, axis=0)
    return js_distance

def main():
    options = get_options()
    infile1 = options.infile1
    infile2 = options.infile2

    df1 = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    df2 = np.loadtxt(infile2, delimiter='\t', dtype='float64')


    js_distance = KDE_JS_divergence(df1, df2)
    print(f"js_distance: {js_distance}")

    #KL_distance = np.sum(rel_entr(z1, z2))
    KL_distance = KDE_KL_divergence(df1, df2)
    print(f"KL_distance: {KL_distance}")

if __name__ == "__main__":
    main()