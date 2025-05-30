import argparse
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

def get_options():
    description = 'Tests importance of classifiers from Pansim gridsearch'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python test_classifiers.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='Directory containing run files from gridsearch run')
    IO.add_argument('--outpref',
                    default = "output",
                    help='Output prefix. Default = "output"')
    IO.add_argument('--params',
                    required=True,
                    help='Comma separated list of predictors to estimate. Must match columns in infile')
    IO.add_argument('--n_estimators',
                type=int,
                default=500,
                help='Number of trees for random forest. Default = 500')

    return parser.parse_args()


def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref
    params = options.params.split(",")
    n_estimators = options.n_estimators
    #print(f"Params-pre: {params}", file=sys.stderr)

    df = pd.read_csv(infile, header=0, sep="\t")
    params = [col for col in params if col in df.columns]
    print(f"Params: {params}", file=sys.stderr)
    print(f"df: {df}", file=sys.stderr)

    # DataFrame with remaining columns
    X = df.drop(columns=params)

    element_dict_list = []
    for element in params:
        print(f"Anaylsing: {element}", file=sys.stderr)
        element_dict = {}
        element_dict["feature"] = element

        y = df[element]
        
        regr = RandomForestRegressor(random_state=42, n_estimators=n_estimators)
        regr.fit(X, y)
        importances = regr.feature_importances_
        
        for i in range(len(importances)):
            element_dict[X.columns[i]] = importances[i]
        
        element_dict_list.append(element_dict)
    
    df_out = pd.DataFrame.from_dict(element_dict_list)
    #print(df_out)
    df_out.to_csv(outpref + "_results.csv", sep = "\t", index = False)


if __name__ == "__main__":
    main()