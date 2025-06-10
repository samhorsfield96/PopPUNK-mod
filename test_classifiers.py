import argparse
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import ElasticNet

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
    IO.add_argument('--regression',
                default="elastic_net",
                choices=['random_forest', 'elastic_net'],
                help='Type of regression to use. Choices: random_forest or random_forest')

    return parser.parse_args()

def clean_dataset(df):
    assert isinstance(df, pd.DataFrame), "df needs to be a pd.DataFrame"
    indices_to_keep = ~df.isin([np.inf, -np.inf]).any(axis=1)
    return df[indices_to_keep].astype(np.float64)

def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref
    params = options.params.split(",")
    n_estimators = options.n_estimators
    regression = options.regression
    #print(f"Params-pre: {params}", file=sys.stderr)

    df = pd.read_csv(infile, header=0, sep="\t").astype(np.float64)
    
    # clean up data
    float32_max = np.finfo(np.float32).max  # ~3.4e+38
    df[df > float32_max] = np.nan  # Replace overflow-prone values with NaN
    if regression == "elastic_net":
        df.dropna(inplace=True)
    
    params = [col for col in params if col in df.columns]
    print(f"Params: {params}", file=sys.stderr)
    print(f"df: {df}", file=sys.stderr)

    # DataFrame with remaining columns
    X = df.drop(columns=params).astype(np.float64)

    print("Max value in X:", np.nanmax(X))
    print("Min value in X:", np.nanmin(X))
    print("Any inf in X:", np.any(np.isinf(X)))
    print("Any NaN in X:", np.any(np.isnan(X)))

    element_dict_list = []
    for element in params:
        print(f"Anaylsing: {element}", file=sys.stderr)
        element_dict = {}
        element_dict["feature"] = element

        y = df[element].astype(np.float64)
        
        if regression == "random_forest":
            regr = RandomForestRegressor(random_state=42, n_estimators=n_estimators)
            regr.fit(X, y)
            importances = regr.feature_importances_

            for i in range(len(importances)):
                element_dict[X.columns[i]] = importances[i]

        else:
            regr = ElasticNet(random_state=0)
            regr.fit(X, y)
            coefficients = regr.coef_
        
            feature_importance = pd.Series(index = X.columns, data = np.abs(regr.coef_))
            total_feature_importance = feature_importance.sum()

            feature_importance = feature_importance / total_feature_importance
            importances = feature_importance.to_dict()
        
            for i in range(len(importances)):
                element_dict[X.columns[i]] = importances[X.columns[i]]
        
        element_dict_list.append(element_dict)
    
    df_out = pd.DataFrame.from_dict(element_dict_list)
    #print(df_out)
    df_out.to_csv(outpref + "_results.csv", sep = "\t", index = False)


if __name__ == "__main__":
    main()