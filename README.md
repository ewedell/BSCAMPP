# BSCAMPP and SCAMPP - Two Scalable Phylogenetic Placement Methods and Frameworks
[![PyPI - Version](https://img.shields.io/pypi/v/bscampp)](https://pypi.org/project/bscampp/#history)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bscampp)](https://pypi.org/project/bscampp/)
[![Build Status](https://img.shields.io/github/actions/workflow/status/ewedell/BSCAMPP/python-package.yml?branch=main&label=build)](https://github.com/ewedell/BSCAMPP/)
[![PyPI - License](https://img.shields.io/pypi/l/bscampp)](https://github.com/ewedell/BSCAMPP/blob/main/LICENSE)
[![Changelog](https://img.shields.io/badge/CHANGELOG-grey)](https://github.com/ewedell/BSCAMPP/blob/main/CHANGELOG.md)

**Table of Contents**
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Example Code and Data](#example-code-and-data)

# Overview
* **Inputs**
  1. Reference tree to place sequences into.
  2. Alignment of reference sequences (protein or nucleotide).
  3. Alignment of query sequences (can be combined with ii.).
  4. Tree info file.
     - (EPA-ng as base method), RAxML-ng info file, typically with suffix `.bestModel`.
     - (pplacer as base method), RAxML-ng or FastTree log file containing model parameters.
* **Output**
  1. Placement results of query sequences in the reference tree in `.jplace` format.


SCAMPP and BSCAMPP are two scalable solutions for phylogenetic placement. SCAMPP is designed more for accuracy
and BSCAMPP is designed more for speed. BSCAMPP achieves some magnitudes of speedup compared to SCAMPP.
The core algorithm is described in detail at <https://doi.org/10.1101/2022.10.26.513936>.
In short, Both frameworks in default use EPA-ng as the base placement method, allowing it to scale to placement trees
of at least ~200,000 leaves. Our two methods achieve this by extracting appropriate subtrees and assigning each query
to its most fitting subtree.

They are divide-and-conquer frameworks and can be used with any base placement methods (e.g., `pplacer` as well).
Currently, BSCAMPP and SCAMPP are implemented with `epa-ng` and `pplacer`.

#### BSCAMPP
It is recommended that BSCAMPP be used with subtrees of size 2000 and with 5 votes based on current best results,
especially if sequences are fragmentary. Defaults for the subtree size and number of votes are set to 2,000 and
5 respectively (see [Usage](#usage) for more details on customizing BSCAMPP).

#### SCAMPP
SCAMPP is also implemented in BSCAMPP, originally from <https://github.com/chry04/PLUSplacer>.
Its default also uses EPA-ng and a subtree size of 2,000.
The user can invoke SCAMPP by running `run_scampp.py` or `scampp` (if installed with PyPI) after installation.

# Installation
BSCAMPP and SCAMPP were tested on **Python 3.8 to 3.12**. There are two ways to install:
(1) with PyPI, or (2) from this GitHub repository. If you have any difficulties installing or running BSCAMPP or SCAMPP,
please contact Eleanor Wedell (ewedell2@illinois.edu).

### External requirements
* **Base placement method**:
  EPA-ng and/or pplacer are requirements since BSCAMPP and SCAMPP will use them as the base phylogenetic placement methods.
  By default, the software will search for binary executables of `pplacer` and `epa-ng` in the user's environment when running for the first time.
  We also included a compiled version of `pplacer` for the Linux system under `bscampp/tools`.
* **(macOS) EPA-ng and pplacer**:
  For macOS systems, you need to have `EPA-ng` and/or `pplacer` installed in your environment.
* **cmake and C++ OpenMP**:
  We also use C++ with OpenMP to speed up the similarity comparison between sequences, which is required to run the pre-compiled binaries.
  In the case of macOS systems, the software will automatically recompile the binaries using `cmake`, which is required to be installed.

### (1) Install with `pip`
The easiest way to install BSCAMPP and SCAMPP is to use `pip install`. This will also install all required Python packages.

```bash
# 1. install with pip (--user if no root access)
pip install bscampp [--user]

# 2. Four binary executables will be installed. The first time
#    running any will create a config file at
#    ~/.bscampp/main.config that resolves the links to all
#    external software (e.g., epa-ng, pplacer)

# ---- BSCAMPP functions
bscampp [-h]    # or
run_bscampp.py [-h]

# ---- SCAMPP functions
scampp  [-h]    # or
run_scampp.py
```

### (2) Install from GitHub
Alternatively, the user can clone this GitHub repository and install the required packages manually.

#### Requirements
```bash
python>=3.7
ConfigParser>=5.0.0
numpy>=1.21.6
treeswift>=1.1.45
taxtastic>=0.9.3
```

```bash
# 1. Close the GitHub repo
git clone https://github.com/ewedell/BSCAMPP.git

# 2. Install all requirements
pip install -r requirements.txt

# 3. Execute BSCAMPP/SCAMPP executables
python run_bscampp.py [-h]
python run_scampp.py [-h]
```

# Usage
All parameter settings can be found by running
```bash
run_bscampp.py -h   #OR
run_scampp.py -h
```

### (1) Default case (`epa-ng`)
```bash
# for BSCAMPP
run_bscampp.py -i [raxml best model] -t [reference tree] -a [alignment file]

# for SCAMPP
run_scampp.py -i [raxml best model] -t [reference tree] -a [alignment file]
```
BSCAMPP and SCAMPP in default mode run EPA-ng as the base method. `[alignment file]` should
contain both sequences from the placement tree and the query sequences to be placed.
This will create an output directory `bscampp_output` and write the placement results to
`bscampp_output/bscampp_result.jplace`.

### (2) Separately giving query alignment and finer control of outputs
```bash
run_bscampp.py -i [raxml best model] -t [reference tree] -a [reference alignment] \
    -q [query sequence alignment] -d [output directory] -o [output name] \
    --threads [num cpus]
```

### (3) Using `pplacer` as the base placement method
```bash
run_bscampp.py -i [logfile from either RAxML/FastTree] -t [reference tree] \
    -a [reference alignment] -q [query sequence alignment] \
    --placement-method pplacer
```
### (4) Changing the number of votes to 15 for BSCAMPP
```bash
run_bscampp.py -i [raxml best model] -t [reference tree] -a [reference alignment] \
    -q [query sequence alignment] -V 15
```

### More comprehensive usage
```bash
> usage: run_bscampp.py [-h] [-v] [--placement-method {epa-ng,pplacer}] -i
>                       INFO_PATH -t TREE_PATH -a ALN_PATH [-q QALN_PATH]
>                       [-d OUTDIR] [-o OUTNAME] [--threads NUM_CPUS] [-m MODEL]
>                       [-b SUBTREESIZE] [-V VOTES]
>                       [--similarityflag SIMILARITYFLAG] [-n TMPFILENBR]
>                       [--fragmentflag FRAGMENTFLAG] [--keeptemp KEEPTEMP]
> 
> This program runs BSCAMPP, a scalable phylogenetic placement framework that scales EPA-ng/pplacer to very large tree placement.
> 
> options:
>   -h, --help            show this help message and exit
>   -v, --version         show program's version number and exit
> 
> BASIC PARAMETERS:
>   These are the basic parameters for BSCAMPP.
> 
>   --placement-method {epa-ng,pplacer}
>                         The base placement method to use. Default: epa-ng
>   -i INFO_PATH, --info INFO_PATH, --info-path INFO_PATH
>                         Path to model parameters. E.g., .bestModel from
>                         RAxML/RAxML-ng
>   -t TREE_PATH, --tree TREE_PATH, --tree-path TREE_PATH
>                         Path to reference tree with estimated branch lengths
>   -a ALN_PATH, --alignment ALN_PATH, --aln-path ALN_PATH
>                         Path for reference sequence alignment in FASTA format.
>                         Optionally with query sequences. Query alignment can
>                         be specified with --qaln-path
>   -q QALN_PATH, --qalignment QALN_PATH, --qaln-path QALN_PATH
>                         Optionally provide path to query sequence alignment in
>                         FASTA format. Default: None
>   -d OUTDIR, --outdir OUTDIR
>                         Directory path for output. Default: bscampp_output/
>   -o OUTNAME, --output OUTNAME
>                         Output file name. Default: bscampp_result.jplace
>   --threads NUM_CPUS, --num-cpus NUM_CPUS
>                         Number of cores for parallelization, default: -1 (all)
>   --cpus-per-job CPUS_PER_JOB
>                         Number of cores to use for each job, default: 2
> 
> ADVANCE PARAMETERS:
>   These parameters control how BSCAMPP and SCAMPP are run. The default values are set based on experiments.
> 
>   -b SUBTREESIZE, --subtreesize SUBTREESIZE
>                         Integer size of the subtree. Default: 2000
>   -V VOTES, --votes VOTES
>                         (BSCAMPP only) Number of votes per query sequence.
>                         Default: 5
>   --similarityflag SIMILARITYFLAG
>                         Boolean, True if maximizing sequence similarity
>                         instead of simple Hamming distance (ignoring gap sites
>                         in the query). Default: True
> 
> MISCELLANEOUS PARAMETERS:
>   -n TMPFILENBR, --tmpfilenbr TMPFILENBR
>                         Temporary file indexing. Default: 0
>   --fragmentflag FRAGMENTFLAG
>                         If queries contains fragments. Default: True
>  --subtreetype SUBTREETYPE
>                         (SCAMPP only) Options for collecting nodes for the
>                         subtree - d for edge weighted distances, n for node
>                         distances, h for Hamming distances. Default: d
>   --keeptemp KEEPTEMP   Boolean, True to keep all temporary files. Default:
                        False
```


# Example Code and Data
Example script and data are provided in this GitHub repository in `examples/`.
The data is originally from the
[RNAsim-VS datasets](https://doi.org/10.1093/sysbio/syz063).
* `examples/run_bscampp.sh`: contains a simple script to test BSCAMPP with
  `epa-ng` or `pplacer`, placing 200 query sequences to a 10000-leaf placement
  tree. The info file is from RAxML-ng when running `epa-ng`, and from
  FastTree-2 when running `pplacer`.
  - `run_bscampp.sh` will invoke BSCAMPP with `epa-ng`.
  - `run_bscampp.sh pplacer` will invoke BSCAMPP with `pplacer`.
* `examples/run_scampp.sh`: the same test script but running SCAMPP.
