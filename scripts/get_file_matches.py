from pathlib import Path
import argparse
import os

def get_options():
    description = 'Matches sample IDs in a large nested directory.'
    
    parser = argparse.ArgumentParser(description=description,
                                     prog='python get_file_matches.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='Tab-separated input file, with species identifier in first column and sample ID in second column')
    IO.add_argument('--outpref',
                    default="output.txt",
                    help='Output prefix. Default = "output"')
    IO.add_argument('--master-dir',
                    required=True,
                    help='Directory of files to search')
    IO.add_argument('--file-ext',
                    required=True,
                    help='File extension to match to e.g. .aln.fas')
    IO.add_argument('--concat',
                    default=False,
                    action="store_true",
                    help='Concatentate contents of found files into single file')
    
    return parser.parse_args()

# search recursively
def find_files_with_extension(root_dir, extension, sample_dict, outpref, concat):
    
    root = Path(root_dir)

    num_searches = len(sample_dict)
    count = 0
    
    found = set()
    with open(outpref + ".tsv", "w") as o1, open(outpref + ".concat", "w") as o2:
        for file in root.rglob(f'*.{extension.lstrip(".")}'):
            base = os.path.basename(file).split(".")[0]
            if base in found:
                continue
            if base in sample_dict:
                o1.write(f"{sample_dict[base]}\t{base}\t{file}\n")
                #print(f"{sample_dict[base]}\t{base}\t{file}\n")

                if concat:
                    with open(file, "r") as i2:
                        o2.write(i2.read())

                found.add(base)
                print(f"Found: {len(found)}")
            
            if len(found) >= num_searches:
                break

    
    #return list(root.rglob(f'*.{extension.lstrip(".")}'))

def main():
    options = get_options()

    infile = options.infile
    outpref = options.outpref
    master_dir = options.master_dir
    concat = options.concat
    file_ext = options.file_ext

    sample_dict = {}
    print("Reading index file...")
    with open(infile, "r") as i1:
        for line in i1:
            taxon_id, sample_id = line.rstrip().split("\t")
            sample_dict[sample_id] = taxon_id

    #print("Generating file list...")
    print("Writing output...")
    find_files_with_extension(master_dir, file_ext, sample_dict, outpref, concat)

    # print("Writing output...")
    # with open(outpref + ".tsv", "w") as o1, open(outpref + ".concat", "w") as o2:
    #     for file in files_list:
    #         base = os.path.basename(file).split(".")[0]
    #         if base in sample_dict:
    #             o1.write(f"{sample_dict[base]}\t{base}\t{file}\n")

    #             if concat:
    #                 with open(file, "r") as i2:
    #                     o2.write(i2.read())

    if not concat:
        if os.path.exists(outpref + ".concat"):
            os.remove(outpref + ".concat")
        else:
            pass

if __name__ == "__main__":
    main()
    