# REcomp2

Pipeline for comparative analysis of potentially unlimited number of RepeatExplorer runs

## Description

Pipeline performs clustering of consensues and 'other' contigs (optionally) from potentially unlimited number of RepeatExplorer runs based on sequences homology that detecting using BLASTn.

Main steps of pipeline's work:  

1. Collection of all sequences for analysis (consensuses, references (if presented) and contigs (optionally)) in one FASTA file

2. Spliting FASTA file in chunks for parallel BLAST

3. All to all chunks BLAST for building of the connectivity table

4. Dropping of junk alignments using K-means or agglomerative clustering

5. Building a graph in memory using QuickUnion algorithm

6. Extraction all connectivity components in graph

7. Filtration of connectivity components based on their content

8. Writing of selected superclusters in separate FASTA

9. Filtering contigs if necessary using BLAST

10. Report generation with information about each selected supercluster

## Installation

Install Anaconda with python 3 from [official website](https://www.anaconda.com/products/individual)

Download and install BLAST+ from [official website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Check the presence of path to installation folder in PATH variable.

```bash
git clone https://github.com/alermol/REcomp2.git
cd REcomp2/
conda env create -f environment.yaml
conda activate recomp
cd REcomp/
./test.py
```

If test run has been done successfully than pipeline is ready to work

Pipeline was tested on Ubuntu 18.04.3

## Usage

The help message and available options can be accessed using

```bash
./REcomp.py -h # short option
./REcomp.py --help # long option
```

which gives the following output

```None
usage: REcomp.py [-h] [-v] [-r REF] [-l] [-c CPU] [-io] [-ir]
                 [--evalue EVALUE] [--low-memory] [-ss {blastn,megablast}]  
                 path prefix out

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
  --low-memory          use small amount of RAM for 'all to all' blast by using small chunk size (1000)
                        but it can take much time (default chunk size: 10000)
  -ss {blastn,megablast}, -superclusters-search {blastn,megablast}
                        alignments for union of sequences in supercluster can be performed either  
                        blastn or megablast (default): blastn is slower and required more RAM
                        but more sensitive
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
If this argument in set `REcomp.log` will be saved in output directory.

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

### `--low-memory`

**Expects**: *None*  
**Default**: *False*  
If this parameters in set, fasta file will be splitted in chunk with less size. This allow to use less amount of memory during 'all to all blast' step.

### `-ss or --superclusters-search`

**Expect**: *STRING (blast or megablast)*  
**Defailt**: *megablast*  
Alignments for union of sequences in superclusters can be performed either `blastn` or `megablast`. `blastn` is slower and required more RAM but more sensitive.

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
