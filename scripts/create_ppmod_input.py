import argparse
import glob
from pathlib import Path

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Creates input file for running of ppmod on real genomes')

    # input options
    parser.add_argument('--indir-pan', help='Path to directory containing output from analyse_pangenome.py',
                                    required=True)
    parser.add_argument('--indir-dist', help='Path to directory containing distances. Should have same naming convention as indir-pan.',
                                    required=True)
    parser.add_argument('--cut-string', help='String to cut from indir-pan entries to get matching name',
                                    default="_core_")
    parser.add_argument('--outpref', help='Output prefix',
                                    default="pp-mod_input.txt")

    return parser.parse_args()

def main():

    # Check input ok
    args = get_options()

    pangenome_files = glob.glob(args.indir_pan + "/*")
    distance_files = glob.glob(args.indir_dist + "/*")

    distance_files_stem = [Path(filename).stem.split(".")[0] for filename in distance_files]

    # match up filenames
    distance_file_order = []
    name_list = []
    for index, filename in enumerate(pangenome_files):
        parsed_filename = Path(filename).stem
        name = parsed_filename.split(args.cut_string)[0]

        name_list.append(name)
        pos = [i for i, s in enumerate(distance_files_stem) if name == s][0]
        distance_file_order.append(pos)
    
    distance_files = [distance_files[i] for i in distance_file_order]

    with open(args.outpref + ".txt", "w") as o:
        o.write("name\tpan_genes\tcore_genes\tintermediate_genes\trare_genes\tavg_gene_freq\tdist_file\n")
        for distfile, panfile, name in zip(distance_files, pangenome_files, name_list):
            pan_dict = {}
            with open(panfile, "r") as file1:
                for line in file1:
                    split_line = line.rstrip().split("\t")
                    pan_dict[split_line[0]] = split_line[1]
            
            o.write(f"{name}\t{pan_dict["pan_genes"]}\t{pan_dict["core_genes"]}\t{pan_dict["intermediate_genes"]}\t{pan_dict["rare_genes"]}\t{pan_dict["avg_gene_freq"]}\t{distfile}\n")


if __name__ == "__main__":
    main()
