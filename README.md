# Reciprocal_introgression
code for simulating and testing reciprocal introgression events

## Intended use of this repo

We will use the repo to:
* Simulate introgression events and generate mulitple sequence alignments from these simulations
* Measure introgression statistics (e.g. D-statistics) on these simulated datasets.


## Running Int_sim.py
* To run Int_sim.py you'll need to install several python modules. We recommend using conda environments to install these needed dependencies
* Example instructions for creating a conda env and installing the needed dependencies. 


```
#Create env
conda create --name int_sim python=3.9.7
conda activate int_sim

#Install modules
conda install -c conda-forge msprime
conda install -c conda-forge matplotlib
conda install -c conda-forge argparse
```

* See the help menu for Int_sim.py to get list of required arguments

```
python Int_sim.py -h
```

