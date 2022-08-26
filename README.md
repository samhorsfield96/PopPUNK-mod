# PopPUNK-mod
PopPUNK modelling: simulates core and accessory genome divergence and calculates Hamming and Jaccard distances, fitting to PopPUNK models.

### Dependencies
- tqdm
- scipy
- numpy
- matplotlib
- scikit-learn
- pandas
- elfi
- numba

### Usage for simulator

```
usage: python run_sim.py [-h] [--core-size CORE_SIZE]
                              [--core-var CORE_VAR]
                              [--base-freq BASE_FREQ]
                              [--base-mu BASE_MU]
                              [--start-gene-freq START_GENE_FREQ]
                              [--avg-gene-freq AVG_GENE_FREQ]
                              [--num-core NUM_CORE] [--num-pan NUM_PAN]
                              [--core-mu CORE_MU] [--acc-mu ACC_MU]
                              [--core-sites CORE_SITES]
                              [--acc-sites ACC_SITES]
                              [--core-gamma-shape CORE_GAMMA_SHAPE]
                              [--core-gamma-scale CORE_GAMMA_SCALE]
                              [--acc-gamma-shape ACC_GAMMA_SHAPE]
                              [--acc-gamma-scale ACC_GAMMA_SCALE]
                              [--core-sites-man CORE_SITES_MAN]
                              [--acc-sites-man ACC_SITES_MAN]
                              [--num-sim NUM_SIM] [--adjust]
                              [--outpref OUTPREF] [--threads THREADS]

Calculate relationship between Hamming/Jaccard distances and core/accessory divergence

optional arguments:
  -h, --help            show this help message and exit

Input/Output options:
  --core-size CORE_SIZE
                        Size of core genome alignment (in bases). Default =
                        1140000
  --core-var CORE_VAR   Number of variant sites in core. Default = 106196
  --base-freq BASE_FREQ
                        Base frequencies in starting core genome in order
                        "A,C,G,T". Default = "0.25,0.25,0.25,0.25"
  --base-mu BASE_MU     Mutation rates from all other bases to each base, in
                        order "A,C,G,T". Default = "0.25,0.25,0.25,0.25"
  --start-gene-freq START_GENE_FREQ
                        Gene frequencies in starting accessory genome in order
                        "0,1". Default = "0.5,0.5"
  --avg-gene-freq AVG_GENE_FREQ
                        Average gene frequency in accessory genome. Default =
                        "0.5"
  --num-core NUM_CORE   Number of core genes. Default = 1194
  --num-pan NUM_PAN     Number of genes in pangenome. Default = 5442
  --core-mu CORE_MU     Range of core genome mutation rate values (mutations
                        per site per genome) in form start,stop,step. Default
                        = "0,2,0.2"
  --acc-mu ACC_MU       Range of accessory gene gain/loss rates (change per
                        gene per genome) in form start,stop,step. Default =
                        "0,2,0.2"
  --core-sites CORE_SITES
                        Number of different core site mutation rates. Default
                        = 3
  --acc-sites ACC_SITES
                        Number of different accessory site mutation rates.
                        Default = 3
  --core-gamma-shape CORE_GAMMA_SHAPE
                        Shape parameter for core per-site substitution rates.
                        Default = 20.0
  --core-gamma-scale CORE_GAMMA_SCALE
                        Scale parameter for core per-site substitution rates.
                        Default = 1.0
  --acc-gamma-shape ACC_GAMMA_SHAPE
                        Shape parameter for accessory per-site substitution
                        rates. Default = 20.0
  --acc-gamma-scale ACC_GAMMA_SCALE
                        Scale parameter for accessory per-site substitution
                        rates. Default = 1.0
  --core-sites-man CORE_SITES_MAN
                        Manual core per-site mutation rates. Must sum to 1.
                        Default = None
  --acc-sites-man ACC_SITES_MAN
                        Manual accessory per-site mutation rates. Must sum to
                        1. Default = None
  --num-sim NUM_SIM     Number of simulations to run. Default = 1
  --adjust              Adjust core and accessory distances for invariant
                        sites. Default = False
  --outpref OUTPREF     Output prefix. Default = "./"
  --threads THREADS     Number of threads. Default = 1
```

### Usage for ELFI model fit

```
usage: python run_ELFI.py [-h] [--core-size CORE_SIZE]
                          [--pan-size PAN_SIZE] 
                          [--max-acc-vs-core MAX_ACC_VS_CORE] 
                          [--num-steps NUM_STEPS] 
                          [--base-mu BASE_MU] 
                          [--avg-gene-freq AVG_GENE_FREQ]
                          [--batch-size BATCH_SIZE] 
                          [--samples SAMPLES] 
                          [--qnt QNT] 
                          [--init-evidence INIT_EVIDENCE] 
                          [--update-int UPDATE_INT] 
                          [--acq-noise-var ACQ_NOISE_VAR] 
                          [--n-evidence N_EVIDENCE]
                          [--data-dir DATA_DIR] 
                          [--data-pref DATA_PREF] 
                          [--seed SEED] 
                          [--summary {quantile,mean}] 
                          [--mode {ABC,BOLFI}] 
                          [--complexity {simple,intermediate}] 
                          [--outpref OUTPREF]
                          [--threads THREADS] 
                          [--cluster]

Fit model to PopPUNK data using Approximate Baysesian computation

options:
  -h, --help            show this help message and exit

Input/Output options:
  --core-size CORE_SIZE
                        Number of positions in core genome. Default = 1000
  --pan-size PAN_SIZE   Number of positions in pangenome. Default = 1000
  --max-acc-vs-core MAX_ACC_VS_CORE
                        Maximum ratio between accessory and core genome evolution. Default = 1000
  --num-steps NUM_STEPS
                        Number of steps to take in increasing divergence. Default = 50
  --base-mu BASE_MU     Mutation rates from all other bases to each base, in order "A,C,G,T". Default = "0.25,0.25,0.25,0.25"
  --avg-gene-freq AVG_GENE_FREQ
                        Average gene frequency in accessory genome. Default = "0.5"
  --batch-size BATCH_SIZE
                        Batch size for processing. Default = 10000
  --samples SAMPLES     No. samples for posterior estimation. Default = 1000
  --qnt QNT             Quantile of the samples with smallest discrepancies is accepted. Default = 0.01
  --init-evidence INIT_EVIDENCE
                        Number of initialization points sampled straight from the priors before starting to optimize the acquisition of points. Default = 5000
  --update-int UPDATE_INT
                        Defines how often the GP hyperparameters are optimized. Default = 10
  --acq-noise-var ACQ_NOISE_VAR
                        Defines the diagonal covariance of noise added to the acquired points. Default = 0.1
  --n-evidence N_EVIDENCE
                        Evidence points requested (including init-evidence). Default = 5000
  --data-dir DATA_DIR   Directory containing popPUNK distance files.
  --data-pref DATA_PREF
                        Prefix of popPUNK distance file(s).
  --seed SEED           Seed for random number generation. Default = 254.
  --summary {quantile,mean}
                        Mode for summary statistics, either "mean" or "quantile". Default = "quantile".
  --mode {ABC,BOLFI}    Mode for running model fit, either "ABC" or "BOLFI". Default = "ABC".
  --complexity {simple,intermediate}
                        Model complexity. If simple, predict only a/pi and gene gain rate. If intermediate, predict prior two and fast gene site mu.Default = "simple".
  --outpref OUTPREF     Output prefix. Default = "./"
  --threads THREADS     Number of threads. Default = 1
  --cluster             Parallelise using ipyparallel if using cluster. Default = False
```