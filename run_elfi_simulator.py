import numpy as np
import argparse
import matplotlib.pyplot as plt
rng = np.random.default_rng()
import subprocess
from run_elfi import read_distfile
import sys

def get_options():
    description = 'Run simulator of gene gain model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_elfi_simulator.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--core_size',
                    type=int,
                    default=1200000,
                    help='Number of positions in core genome. Default = 1200000 ')
    IO.add_argument('--pan_size',
                    type=int,
                    default=6000,
                    help='Number of positions in pangenome. Default = 6000 ')
    IO.add_argument('--core_mu',
                    type=float,
                    default=0.05,
                    help='Maximum pairwise distance for core genome. Default = 0.05 ')
    IO.add_argument('--pan_mu',
                    type=float,
                    default=0.05,
                    help='Maximum pairwise distance for accessory genome. Default = 0.05 ')
    IO.add_argument('--proportion_fast',
                    type=float,
                    default=0.5,
                    help='Proportion of pangenome that is evolves quickly. Default = 0.5')
    IO.add_argument('--speed_fast',
                    type=float,
                    default=2.0,
                    help='Fold-speed increase of pangenome that is evolves quickly over that which evolves slowly. Default = 2.0')
    IO.add_argument('--pop_size',
                    type=int,
                    default=1000,
                    help='Population size for Wright-Fisher model. Default = 1000 ')
    IO.add_argument('--n_gen',
                    type=int,
                    default=100,
                    help='Number of generations for Wright-Fisher model. Default = 100 ')
    IO.add_argument('--avg_gene_freq',
                    type=float,
                    default=0.5,
                    help='Average gene frequency in accessory genome. '
                         'Default = "0.5" ')
    IO.add_argument('--max_distances',
                    type=int,
                    default=100000,
                    help='Number of distances to sample with Pansim. Default = 100000')
    IO.add_argument('--outpref',
                    default="sim",
                    help='Output prefix. Default = "sim"')
    IO.add_argument('--seed',
                    type=int,
                    default=254,
                    help='Seed for random number generation. Default = 254. ')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1')
    IO.add_argument('--pansim_exe',
                    required=True,
                    help='Path to pansim executable.')

    return parser.parse_args()

if __name__ == "__main__":
    # core_mu = 0.05
    # pan_mu = 0.05
    # proportion_fast = 0.5
    # speed_fast = 2.0
    # core_size = 1200000
    # pan_size = 6000
    # avg_gene_freq = 0.5
    # n_gen = 10
    # pop_size = 1000
    # max_distances = 100000
    # threads = 8
    # outpref = "test"
    # pansim_exe = "/home/shorsfield/software/Pansim/pansim/target/release/pansim"
    # seed = 254

    options = get_options()
    core_mu = options.core_mu
    pan_mu = options.pan_mu
    proportion_fast = options.proportion_fast
    speed_fast = options.speed_fast
    core_size = options.core_size
    pan_size = options.pan_size
    avg_gene_freq = options.avg_gene_freq
    n_gen = options.n_gen
    pop_size = options.pop_size
    max_distances = options.max_distances
    threads = options.threads
    outpref = options.outpref
    pansim_exe = options.pansim_exe
    seed = options.seed

    command = pansim_exe + ' --avg_gene_freq {avg_gene_freq} --pan_mu {pan_mu} --proportion_fast {proportion_fast} --speed_fast {speed_fast} --core_mu {core_mu} --seed {seed} --pop_size {pop_size} --core_size {core_size} --pan_size {pan_size} --n_gen {n_gen} --max_distances {max_distances} --output {output_filename} --threads {threads}'.format(avg_gene_freq=avg_gene_freq, pan_mu=pan_mu, proportion_fast=proportion_fast, speed_fast=speed_fast, core_mu=core_mu, seed=seed, pop_size=pop_size, core_size=core_size, pan_size=pan_size, n_gen=n_gen, max_distances=max_distances, output_filename=outpref + ".txt", threads=threads)

    print("Simulating...")
    try:
        result = subprocess.run(command.split(" "))
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print(e.output)
        sys.exit(1)
    
    print("Generating figures...")
    df = read_distfile(outpref + ".txt")

    fig, ax = plt.subplots()

    x = df["Core"]
    y = df["Accessory"]

    ax.scatter(x, y, s=10, alpha=0.3)

    ax.set_xlabel("Core distance")
    ax.set_ylabel("Accessory distance")

    fig.savefig(outpref + "_sim" + ".png")


