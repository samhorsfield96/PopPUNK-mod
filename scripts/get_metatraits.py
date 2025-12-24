import pandas as pd
import requests
import time
import tqdm

def get_options():
    description = 'Download trait data from metatraits (https://metatraits.embl.de/).'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python get_metatraits.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--metadata',
                    required=True,
                    help='Input metadata e.g. https://metatraits.embl.de/static/downloads/GTDB2NCBI.tsv.gz.')
    IO.add_argument('--outpref',
                    default="output",
                    help='Output prefix.')
    IO.add_argument('--database',
                    default="habitat_clustering:SPIREv1",
                    choices=['BacDive:2025-07-29', 'BacDive-AI:proGenomes3', 'BacDive-AI:SPIREv1', 'BV-BRC:2025-07-29', 'GenomeSPOT:proGenomes3', 'GenomeSPOT:SPIREv1', 'Genome statistics:proGenomes3', 'Genome statistics:SPIREv1', 'habitat_clustering:SPIREv1', 'JGI/GOLD:2025-07-29', 'JGI/GOLD,JGI/IMG:2025-07-29', 'MICROPHERRET:proGenomes3', 'MICROPHERRET:SPIREv1', 'Traitar:proGenomes3', 'Traitar:SPIREv1'],
                    help='Database to pull from. Default = "habitat_clustering:SPIREv1".')
    IO.add_argument('--base_url',
                default="https://metatraits.embl.de/api/v1/traits/taxonomy",
                help='Base URL for search. Default = ""https://metatraits.embl.de/api/v1/traits/taxonomy"')
    
    return parser.parse_args()

def main():
    options = get_options()

    # Load metadata
    metadata_file = options.metadata
    outpref = options.outpref
    databases_param = options.database
    base_url = option.base_url

    df = pd.read_csv(metadata_file, sep='\t')

    # Base API URL
    base_url += "/{}"

    results = []

    for index, row in tqdm.tqdm(df.iterrows(), total=len(df), desc="Querying API"):
        ncbi_taxon_id = row['taxonID NCBI']
        
        # Build URL with database filter
        url = base_url.format(ncbi_taxon_id)
        params = {
            'databases': databases_param
        }
        
        try:
            response = requests.get(url, params=params)
            if response.status_code == 200:
                traits_json = response.json()
                
                trait_dict = {
                    'taxonID NCBI': ncbi_taxon_id,
                    'species (NCBI)': row.get('species (NCBI)', ''),
                }
                
                for feature in traits_json:
                    feat_name = feature.get('feature')
                    feat_value = feature.get('value')
                    trait_dict[feat_name] = feat_value
                
                results.append(trait_dict)
            else:
                print(f"Warning: Status {response.status_code} for taxon ID {ncbi_taxon_id}")
                results.append({'taxonID NCBI': ncbi_taxon_id, 'species (NCBI)': row.get('species (NCBI)', '')})
        except Exception as e:
            print(f"Error querying taxon ID {ncbi_taxon_id}: {e}")
            results.append({'taxonID NCBI': ncbi_taxon_id, 'species (NCBI)': row.get('species (NCBI)', '')})
        
        time.sleep(0.01)


    results_df = pd.DataFrame(results)

    output_df = pd.merge(df, results_df, on=['taxonID NCBI', 'species (NCBI)'], how='left')

    output_df.to_csv(outpref + ".tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()
