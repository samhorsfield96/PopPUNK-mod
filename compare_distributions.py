import argparse
from run_elfi import read_distfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
import seaborn as sns
import pandas as pd

def get_options():
    description = 'Run simulator of gene gain model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_elfi_simulator.py')

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
            help='Output prefix. Default = "joint')

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

    df1 = read_distfile(file1)
    df2 = read_distfile(file2)

    max_acc = max(float(df1["Accessory"].max()), float(df2["Accessory"].max()))
    min_acc = max(float(df1["Accessory"].min()), float(df2["Accessory"].min()))

    num_bins = 1000
    #hist_range = (min_acc, max_acc)
    hist_range = (0, 1)

    # first distribution is x, second is y
    x = np.histogram(df1['Accessory'].to_numpy(), bins=num_bins, range=hist_range)
    y = np.histogram(df2['Accessory'].to_numpy(), bins=num_bins, range=hist_range)

    js_dist = distance.jensenshannon(x[0], y[0])
    print("JS distance: {}".format(str(js_dist)))

    plt.hist(df1["Accessory"], num_bins, alpha=0.5, label='original')
    plt.hist(df2["Accessory"], num_bins, alpha=0.5, label='new')
    plt.legend(loc='upper right')

    # Save the plot
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