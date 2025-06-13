
import numpy as np
import argparse
from scipy.spatial.distance import jensenshannon, canberra, cosine, euclidean
from scipy.special import rel_entr
from run_elfi import get_grid
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
def get_grid_nonsymmetric(min_x, max_x, min_y, max_y, resolution):
    x = np.linspace(min_x, max_x, resolution)
    y = np.linspace(min_y, max_y, resolution)
    xx, yy = np.meshgrid(x, y)
    xy = np.vstack([xx.ravel(), yy.ravel()]).T

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

def create_KDE_dist(df, grid_params):
    scale = np.amax(df, axis = 0)
    df /= scale

    kde = get_kde(df)
    z = generate_samples(grid_params, kde)
    #print(f"z: {z}")
    z = z.ravel()
    z /= z.sum()

    return z

def main():
    options = get_options()
    infile1 = options.infile1
    infile2 = options.infile2

    infile1_df = np.loadtxt(infile1, delimiter='\t', dtype='float64')
    infile2_df = np.loadtxt(infile2, delimiter='\t', dtype='float64')

    scaler = MinMaxScaler()
    scaler.fit(np.vstack([infile1_df, infile2_df]))
    #print(f"infile1_df_pre: {infile1_df}")
    infile1_df = scaler.transform(infile1_df)
    infile2_df = scaler.transform(infile2_df)

    #print(f"infile1_df_post: {infile1_df}")

    max_infile1 = np.max(infile1_df, axis=0)
    min_infile1 = np.min(infile1_df, axis=0)
    max_infile2 = np.max(infile2_df, axis=0)
    min_infile2 = np.min(infile2_df, axis=0)

    # print(f"max_infile1: {max_infile1}")

    #print(f"max_acc1: {max_infile1[1]}, max_core1: {max_infile1[0]}")
    #print(f"max_acc2: {max_infile2[1]}, max_core2: {max_infile2[0]}")
    #print(f"min_acc1: {min_infile1[1]}, min_core1: {min_infile1[0]}")
    #print(f"min_acc2: {min_infile2[1]}, min_core2: {min_infile2[0]}")

    min_x = min(min_infile1[0], min_infile2[0])
    max_x = max(max_infile1[0], max_infile2[0])
    min_y = min(min_infile1[1], min_infile2[1])
    max_y = max(max_infile1[1], max_infile2[1])

    # print(f"min_x: {min_x}")
    # print(f"max_x: {max_x}")
    # print(f"min_y: {min_y}")
    # print(f"max_y: {max_y}")

    #grid_params = get_grid_nonsymmetric(min_x, max_x, min_y, max_y, 100)
    grid_params = get_grid(0, 1, 100)

    #print(f"grid_params[0]: {grid_params[0]}")
    #print(f"grid_params[1]: {grid_params[1]}")
    #print(f"grid_params[2]: {grid_params[2]}")

    #print("z1")
    z1 = create_KDE_dist(infile1_df, grid_params)
    #print("z2")
    z2 = create_KDE_dist(infile2_df, grid_params)

    js_distance = jensenshannon(z1, z2, axis=0)
    print(f"js_distance: {js_distance}")

    #js_distance2 = js_distance_integral_kde(get_kde(infile1_df), get_kde(infile2_df))
    #print(f"js_distance2: {js_distance2}")

    # infile1_arr = np.array([0, min_infile1[0], min_infile1[1], max_infile1[0], max_infile1[1]])
    # infile2_arr = np.array([js_distance, min_infile1[0], min_infile1[1], max_infile1[0], max_infile1[1]])

    # canberra_distance = canberra(infile1_arr, infile2_arr)
    # print(f"canberra_distance: {canberra_distance}")
    # cosine_distance = cosine(infile1_arr, infile2_arr)
    # print(f"cosine_distance: {cosine_distance}")
    # euclidean_distance = euclidean(infile1_arr, infile2_arr)
    # print(f"euclidean_distance: {euclidean_distance}")

if __name__ == "__main__":
    main()