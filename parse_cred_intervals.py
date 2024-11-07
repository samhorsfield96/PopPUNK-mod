import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

def get_options():
    description = 'Run simulator of gene gain model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python run_elfi_simulator.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                required=True,
                help='Input file.')
    IO.add_argument('--name',
                    default=None,
                    help='Name for parameter set. Default is input file name.')
    IO.add_argument('--outfile',
                    default="cred_intervals.txt",
                    help='Output file. Default = "cred_intervals.txt"')

    return parser.parse_args()

def main():
    options = get_options()

    infile = options.infile
    outfile = options.outfile
    name = options.name

    skip = True
    param_dict = {}
    with open(infile, "r") as f1:
        while True:
            line = f1.readline()
            if not line:
                break
            
            split_line = line.split()

            # if at parameter line, then read following lines
            if len(split_line) > 0:
                if split_line[0] == "Parameter":
                    skip = False
                    continue

            if len(split_line) > 0 and skip == False:
                # get parameter name
                parameter = split_line[0][:-1]
                mean = float(split_line[1])
                cred_1 = float(split_line[2])
                cred_2 = float(split_line[3])

                param_dict[parameter] = (mean, cred_1, cred_2)
    
    if name == None:
        name = os.path.splitext(os.path.basename(infile))[0]

    if os.path.isfile(outfile):
        with open(outfile, "a") as o:
            for param, param_val in param_dict.items():
                o.write(str(name) + "\t" + str(param) + "\t" + str(param_val[0]) + "\t" + str(param_val[1]) + "\t" + str(param_val[2]) + "\n")
    else:
        with open(outfile, "w") as o:
            o.write("Name\tParam\tMean\tCred_2.5\tCred_97.5\n")
            for param, param_val in param_dict.items():
                o.write(str(name) + "\t" + str(param) + "\t" + str(param_val[0]) + "\t" + str(param_val[1]) + "\t" + str(param_val[2]) + "\n")


if __name__ == "__main__":
    main()