import argparse
from pathlib import Path

def get_options():
    description = 'Merges gridsearch runs into single file'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python merge_gridsearch_params.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--indir',
                    required=True,
                    help='Directory containing run files from gridsearch run.')
    IO.add_argument('--params',
                    required=True,
                    help='Parameter file for gridsearch run.')
    IO.add_argument('--outpref',
                    default = "output",
                    help='Output prefix. Default = "output"')         

    return parser.parse_args()


def main():
    options = get_options()
    indir = options.indir
    params = options.params
    outpref = options.outpref

    params_list = []
    header_params = None

    with open(params, "r") as i1:
        header_params = i1.readline().rstrip()
        for line in i1:
            params_list.append(line.rstrip())

    #print(params_list[0])
    #print(params_list[-1])

    files = []
    for path in Path(indir).glob("**/*.tsv"):
        # Print the path (file or directory) to the console
        files.append(str(path))

    #print(files)

    header_results = None
    results_dict = {}
    for file in files:
        with open(file, "r") as i1:
            header_results = i1.readline().rstrip()
            for line in i1:
                split_line = line.rstrip().split("\t")
                filepref = split_line[0]

                results_dict[filepref] = line.rstrip()

    #print(results_dict[1000])

    with open(outpref + ".tsv", "w") as o1:
        o1.write(header_results + "\t" + header_params + "\n")

        for filepref, entry in results_dict.items():
            file_idx = int(filepref.split("_")[-1])
            #print(file_idx)
            o1.write(entry + "\t" + params_list[file_idx - 1] + "\n")
            

if __name__ == "__main__":
    main()