# distance_sim
Simulates core and accessory genome divergence and calculates Hamming and Jaccard distances.

### Dependencies
- tqdm
- scipy
- numpy
- matplotlib
- scikit-learn

### Usage

```
usage: python __main__.py [-h] [--core-size CORE_SIZE] [--core-var CORE_VAR] [--base-freq BASE_FREQ] [--base-mu BASE_MU] 
                               [--start-gene-freq START_GENE_FREQ] [--avg-gene-freq AVG_GENE_FREQ] [--num-core NUM_CORE]
                               [--num-pan NUM_PAN] [--core-mu CORE_MU] [--acc-mu ACC_MU] [--num-sim NUM_SIM] [--adjust] 
                               [--outpref OUTPREF] [--threads THREADS]

Calculate relationship between Hamming/Jaccard distances and core/accessory divergence

options:
  -h, --help            show this help message and exit

Input/Output options:
  --core-size CORE_SIZE
                        Size of core genome alignment (in bases). Default = 1140000
  --core-var CORE_VAR   Number of variant sites in core. Default = 106196
  --base-freq BASE_FREQ
                        Base frequencies in starting core genome in order "A,C,G,T". Default = "0.25,0.25,0.25,0.25"
  --base-mu BASE_MU     Mutation rates from all other bases to each base, in order "A,C,G,T". Default = "0.25,0.25,0.25,0.25"
  --start-gene-freq START_GENE_FREQ
                        Gene frequencies in starting accessory genome in order "0,1". Default = "0.5,0.5"
  --avg-gene-freq AVG_GENE_FREQ
                        Average gene frequency in accessory genome. Default = "0.5"
  --num-core NUM_CORE   Number of core genes. Default = 1194
  --num-pan NUM_PAN     Number of genes in pangenome. Default = 5442
  --core-mu CORE_MU     Range of core genome mutation rate values (mutations per site per genome) in form start,stop,step. Default = "0,2,0.2"
  --acc-mu ACC_MU       Range of accessory gene gain/loss rates (change per gene per genome) in form start,stop,step. Default = "0,2,0.2"
  --num-sim NUM_SIM     Number of simulations to run. Default = 1
  --adjust              Adjust core and accessory distances for invariant sites. Default = False
  --outpref OUTPREF     Output prefix. Default = "./"
  --threads THREADS     Number of threads. Default = 1
```