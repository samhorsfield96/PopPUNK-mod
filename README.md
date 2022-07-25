# distance_sim
Simulates core and accessory genome divergence and calculates Hamming and Jaccard distances.

### Dependencies
- tqdm
- scipy
- numpy
- matplotlib

### Usage

```
usage: python __main__.py [-h] [--aln-core ALN_CORE] [--num-core NUM_CORE] [--num-pan NUM_PAN] [--core-mu CORE_MU]
                          [--acc-mu ACC_MU] [--outpref OUTPREF] [--threads THREADS]

Calculate relationship between Hamming/Jaccard distances and core/accessory divergence

options:
  -h, --help           show this help message and exit

Input/Output options:
  --aln-core ALN_CORE  Size of core genome alignment (in bases). Default = 1140000
  --num-core NUM_CORE  Number of core genes. Default = 1194
  --num-pan NUM_PAN    Number of genes in pangenome. Default = 5442
  --core-mu CORE_MU    Range of core genome mutation rate values (mutations per site per genome) in form
                       start,stop,step. Default = "0,2,0.2"
  --acc-mu ACC_MU      Range of accessory gene gain/loss rates (change per gene per genome) in form start,stop,step.
                       Default = "0,2,0.2"
  --outpref OUTPREF    Output prefix. Default = "./"
  --threads THREADS    Number of threads. Default = 1
```