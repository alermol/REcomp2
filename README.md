# REcomp2

Pipeline for comparative analysis potentially unlimited results of RepeatExplorer runs

## Description

## Installation

Pipeline was tested on Ubuntu 18.04.3

Install Anaconda with python 3 from [official website](https://www.anaconda.com/products/individual)

Download and install BLAST+ from [official website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Check the presence of path to installation folder in PATH variable.

Download ZIP-archive

Unzip archive in installation directory

```bash
cd REcomp2-master
conda env create -f environment.yaml
conda activate recomp
pip install python_algorithms
cd REcomp
./test.py
```

If test run has been done successfully than pipeline is ready to work

## Usage

The help message and available options can be accessed using

```bash
./REcomp.py -h # short option
./REcomp.py --help # long option
```

which gives the following output

```None
usage: REcomp.py [-h] [-v] [-r REF] [-l] [-c CPU] [-io] [-ir]
                 [--evalue EVALUE] [--ident IDENTITY_PERC]
                 [--qcov QUERY_COVER] [--low-memory] path prefix out

positional arguments:
  path                  path(s) to RE results (top level)
  prefix                prefix(es) for each paths
  out                   path to output directory

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show version
  -r REF, --references REF
                        path to fasta with references repeats
  -l                    save logfile in output directory
  -c CPU                number of CPU to use
  -io, --include-other  include 'other' contigs and clusters in analysis (default: False)
  -ir, --include-ribosomal
                        include rDNA clusters (rank 4) in analysis (default: False)
  --evalue EVALUE       evalue threshold for alignments for supercluster assembly (default: 1e-05)
  --ident IDENTITY_PERC
                        identity percent threshold for alignment for superclusters assembly (default: 90.0)
  --qcov QUERY_COVER    query cover threshold for alignment for superclusters assembly (default: 80.0)
  --low-memory          use small amount of RAM for 'all to all' blast by using small chunk size (1000)
                        but it can take much time (default chunk size: 10000)
```

The details of each option are given below:

### `path`

**Expects**: *STRING*  
**Default**: *None*  
The input for this argument is string in quotes that contains space separated paths to top level of RepeatExplorer results. Paths must not be repeated. This argument is required.  
**Example**  
```"~/RE_result1 ~/RE_result2"```

### `prefix`

**Expects**: *STRING*  
**Default**: *None*  
This argument take in a string in quotes that contains space seaparted prefixes for each paths in `path` argument. Each one prefix must be matching only one path. Prefixes must not be repeated. Short prefixes are recomended. This argument is requered.  
**Example**  
```"p1 p2"```

### `out`

**Expect**: *STRING (to be used as a path to directory)*  
**Default**: *None*  
This is the last required argument for REcomp2. Points to directory for results saving. If directory does not exist, it will be created.  

### `-r or --references`

**Expects**: *STRING (to be used as a path to file)*  
**Default**: *None*  
File containing references in FASTA format. Files with multiple sequences are also valid. Multi-FASTA must not contain identical records id.

### `-l`

**Expects**: *None*  
**Default**: *False*  
If this argument in set `logfile.log` will be saved in output directory.

### `-c`

**Expects**: *INTEGER*  
**Default**: *All system cores*  
Number of cores for REcomp2 work.

### `-io or --include-other`

**Expects**: *None*  
**Default**: *False*  
By default the 'other' clusters and contigs does not included in analysis. If this arguments in set, contigs from 'other' clusters is included in analysis that can dramatically increase worktime, especially if `--low-memory` (see below) is turn on.

### `-ir or --include-ribosomal`

**Expects**: *None*  
**Default**: *False*  
This parameter activates including in analysis the rDNA consensuses (rank 4).

### `--evalue`

**Expects**: *FLOAT >= 0.0*  
**Default**: *1e-05*  
E-value threshold for BLAST alignment. Alignments with E-value greater than this threshold are not counted.

### `--ident`

**Expects**: *0.0 <= FLOAT <= 100.0*  
**Default**: *90.0*  
Identity percent thershold for BLAST alignment. Alignments with identity percent less than thershold are not counted.

### `--qcov`

**Expects**: *0.0 <= FLOAT <= 100.0*  
**Default**: *80.0*  
Query cover percent thershold for BLAST alignment. Alignments with query cover less than thershold are not counted.

### `--low-memory`

**Expects**: *None*  
**Default**: *False*  
If this parameters in set, fasta file will be splitted in chunk with less size. This allow to use less amount of memory during 'all to all blast' step.

### `-v or --version`

Prints the version info of REcomp2

## Examples


```bash
# find all superclusters using default parameters
./REcomp.py '~/RE_result1 ~/RE_result2' 'p1 p2' ~/REcopm2_output -r ~/references.fasta
# find all superclusters using default parameters and log-file
./REcomp.py '~/RE_result1 ~/RE_result2' 'p1 p2' ~/REcopm2_output -r ~/references.fasta -l
# find all superclusters incliding 'other' using default parameters 
./REcomp.py '~/RE_result1 ~/RE_result2' 'p1 p2' ~/REcopm2_output -r ~/references.fasta -io
# find all superclusters incliding 'other' using less memory
./REcomp.py '~/RE_result1 ~/RE_result2' 'p1 p2' ~/REcopm2_output -r ~/references.fasta -io --low-memory
```
