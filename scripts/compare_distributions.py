import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.optimize import curve_fit

import seaborn as sns
import pandas as pd
from run_elfi import wasserstein_distance, js_distance, rmse, asymptotic_curve

def read_distfile(filename):
    # read first line, determine if csv
    with open(filename, "r") as f:
        first_line = f.readline()
        if "," in first_line:
            obs = pd.read_csv(filename, index_col=None, header=None, sep=",")
        else:
            obs = pd.read_csv(filename, index_col=None, header=None, sep="\t")

    if len(obs.columns) == 2:
        obs.rename(columns={obs.columns[0]: "Core",
                           obs.columns[1]: "Accessory"}, inplace=True)
    elif len(obs.columns) == 4:
        # rename columns
        obs.rename(columns={obs.columns[0]: "Sample1", obs.columns[1] : "Sample2", obs.columns[2]: "Core",
                           obs.columns[3]: "Accessory"}, inplace=True)
    else:
        print("Incorrect number of columns in distfile. Should be 2 or 4.")
        sys.exit(1)

    obs['Core'] = pd.to_numeric(obs['Core'])
    obs['Accessory'] = pd.to_numeric(obs['Accessory'])

    return obs

def get_options():
    description = 'Calculates an assortment of statistics between two 2-dimensional distributions.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python compare_distributions.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--file1',
                    type=str,
                    required=True,
                    help='Input file 1.')
    IO.add_argument('--file2',
                type=str,
                required=True,
                help='Input file 2.')
    IO.add_argument('--outpref',
            type=str,
            default="joint",
            help='Output prefix. Default = "joint"')
    IO.add_argument('--plot',
            action="store_true",
            default=False,
            help='Generate plots')

    return parser.parse_args()

def main():
    #file1 = "/home/shorsfield/software/PopPUNK-mod/simulations/sim_pm_0.05_pf_0.5_sf_2.0.txt"
    #file2 = "/home/shorsfield/software/PopPUNK-mod/simulations/sim_pm_0.025_pf_0.5_sf_2.0.txt"
    #outpref = "/home/shorsfield/software/PopPUNK-mod/sim_pm_0.05_pf_0.5_sf_2.0"
    #outpref = "joint"

    options = get_options()
    file1 = options.file1
    file2 = options.file2
    outpref = options.outpref
    plot = options.plot

    df1 = read_distfile(file1)
    df2 = read_distfile(file2)

    # calculate RMSE between distribution fits
    # fit asymptotic curve
    popt, pcov = curve_fit(asymptotic_curve, df1["Core"], df1["Accessory"], p0=[1.0, 1.0, 0.0])
    
    # pull out values for distribution
    a, b, c = popt

    acc_fit1 = asymptotic_curve(df1['Core'], a, b, c)
    RMSE1 = rmse(df1['Accessory'], acc_fit1)

    acc_fit2 = asymptotic_curve(df2['Core'], a, b, c)
    RMSE2 = rmse(df2['Accessory'], acc_fit2)

    print("RMSE original: {}".format(str(RMSE1)))
    print("RMSE new: {}".format(str(RMSE2)))
    print("RMSE diff: {}".format(str(RMSE2 - RMSE1)))

    num_bins = 500
    #hist_range = (min_acc, max_acc)
    #hist_range = (0, 1)

    #max_core = max(df1['Core'].max(), df2['Core'].max())
    #max_acc = max(df1['Accessory'].max(), df2['Accessory'].max())
    max_core = df1['Core'].max()
    max_acc = df1['Accessory'].max()

    # process distributions
    core_1 = np.histogram(df1['Core'].to_numpy(), bins=num_bins, range=(0, max_core))[0]
    acc_1 = np.histogram(df1['Accessory'].to_numpy(), bins=num_bins, range=(0, max_acc))[0]
    dist_1 = (core_1, acc_1)

    core_2 = np.histogram(df2['Core'].to_numpy(), bins=num_bins, range=(0, max_core))[0]
    acc_2 = np.histogram(df2['Accessory'].to_numpy(), bins=num_bins, range=(0, max_acc))[0]
    dist_2 = (core_2, acc_2)

    print("Calculating distances...")
    wass_dist = wasserstein_distance(dist_1, dist_2)
    print("Wasserstein distance: {}".format(str(wass_dist)))

    js_core = js_distance(dist_1, 0, dist_2)
    js_acc = js_distance(dist_1, 1, dist_2)
    js_avg = (js_core + js_acc) / 2
    print("JS distance core: {}".format(str(js_core)))
    print("JS distance acc: {}".format(str(js_acc)))
    print("JS distance avg: {}".format(str(js_avg)))

    plt.hist(df1["Accessory"], num_bins, alpha=0.5, label='original')
    plt.hist(df2["Accessory"], num_bins, alpha=0.5, label='new')
    plt.legend(loc='upper right')

    # Save the plot
    if plot:
        try:
            plt.savefig(outpref + '_hist.png', dpi=300, bbox_inches='tight')
            print(f"Plot saved as {outpref}_hist_.png")
        except Exception as e:
            print("Error saving the plot:", e)

        plt.close()
        
        joint_df = pd.DataFrame({'x': df1["Accessory"], 'y': df2["Accessory"]})
        #joint_df = pd.DataFrame({'x': x, 'y': y})
        sns.jointplot(data=joint_df, x="x", y="y", alpha = 0.3)

        try:
            plt.savefig(outpref + '_jointplot.png', dpi=300, bbox_inches='tight')
            print(f"Plot saved as {outpref}_jointplot.png")
        except Exception as e:
            print("Error saving the plot:", e)


if __name__ == "__main__":
    main()