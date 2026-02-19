# opsintools

This project contains a set of tools for analysis of opsin structures and sequences.

Each tool has an API and a CLI. Keep in mind that the interfaces might change in the future since the project is under active development.

## Try it

Some of the tools can be run from the browser via [Google Colab](https://colab.research.google.com/):

* [opsinmaphmm](https://colab.research.google.com/github/BejaLab/opsintools/blob/main/colab/opsinmaphmm.ipynb)
* [opsinpdb](https://colab.research.google.com/github/BejaLab/opsintools/blob/main/colab/opsinpdb.ipynb)

## Install it

Installation is currently available for Unix-like systems. There are python and non-python dependencies: the most straightforward installation is via conda which takes care of both.

### conda

Use the provided environment definition: `conda env create -n opsintools -f env.yaml --channel-priority flexible` and then `conda activate opsintools`. Use `mamba` if available.

### pypi

`opsintools` and its python dependencies can be installed with `pip install git+https://github.com/BejaLab/opsintools`. Non-python dependencies should then be installed separately. As of now these include:

* [`dssp`](https://swift.cmbi.umcn.nl/gv/dssp/index.html) (v.3.0.0 tested)
* [`t-coffee`](https://tcoffee.crg.eu/) (v.13.46 tested)
* [`US-align`](https://github.com/pylelab/USalign) (v.20241201 tested)

## The tools

### opsinmap3d

```
usage: opsinmap3d [-h] -i INPUT -d DATADIR -o OUTPUT [-f] [-n N_REPS] [-t THREADS] [--only-exptl] [--prefer-exptl] [--pad-n PAD_N] [--pad-c PAD_C]
                  [--max-seq-id MAX_SEQ_ID] [--methods METHODS]

Opsin homology based on 3D structure.

options:
  -h, --help            show this help message and exit

Arguments:
  -i INPUT              query PDB structure (only chain A will be used)
  -d DATADIR            opsin data directory
  -o OUTPUT             output directory
  -f                    whether to overwrite files in the output directory if it exists
  -n N_REPS             maximum number of template structures to use (default: 20)
  -t THREADS            number of threads to use (default: 1)

Advanced:
  --only-exptl          only use experimental structures (default: use predictions as well)
  --prefer-exptl        prefer experimental structure over predictions (default: choose by similarity only)
  --pad-n PAD_N         N-terminal padding for query trimming (default: 30)
  --pad-c PAD_C         C-terminal padding for query trimming (default: 30)
  --max-seq-id MAX_SEQ_ID
                        maximum sequence identity for predicted templates (default: 0.9)
  --methods METHODS     t-coffee methods to use (default: sap_pair,mustang_pair,t_coffee_msa,probcons_msa)
```

Opsin homology based on protein structures. The input is an opsin structure in PDB format and a suitable reference dataset (can be dowloaded from [opsintools-build](https://github.com/BejaLab/opsintools-build/releases)), the output is the mapping between positions in the query and the positions in the reference in the transmembrane regions. Both the CLI and the API generate output files located in the specified directory:

* `aln_to_ref.txt` - alignment of the query to the reference
* `trimmed.pdb` - trimmed query structure
* `t_coffee.aln` - structural alignment
* `t_coffee.aln.score_ascii` - alignment scores
* `t_coffee.aln.log` - structural alignment log
* `opsinmap.json` - json file with the position mapping

The API function `opsinmap3d` returns the dictionary mapping reference positions to the query positions as its output.

See `opsinmap3d -h` or `from opsintools import opsinmap3d; help(opsinmap3d)` for more details.

### opsinalign3d

```
usage: opsinalign3d [-h] -i INPUT [INPUT ...] -o OUTPUT [-f] [-t THREADS] [--methods METHODS]

Opsin alignment with t-coffee.

options:
  -h, --help            show this help message and exit

Arguments:
  -i INPUT [INPUT ...]  input PDB structures (only chain A will be used)
  -o OUTPUT             output directory
  -f                    whether to overwrite files in the output directory if it exists
  -t THREADS            number of threads to use (default: 1)

Advanced:
  --methods METHODS     t-coffee methods to use (default: sap_pair,mustang_pair,t_coffee_msa,probcons_msa)
```

Runs `t-coffee` on a set of PDB files and parses the results. The output is an object of the custom class `Tcoffee`. The output directory will contain:

* `t_coffee.aln` - structural alignment
* `t_coffee.log` - structural alignment log

See `opsinalign3d -h` or `from opsintools import opsinalign3d; help(opsinalign3d)` for more details.

Currently supported methods are: `sap_pair`, `mustang_pair`, `t_coffee_msa`, `probcons_msa`, `vmaffteinsi_msa`, `mafftfftns1_msa`, `mafftfftnsi_msa`, `mafftginsi_msa`, `mafftlinsi_msa`, `mafft_msa`, `mafftnwnsi_msa`, `mafftsparsecore_msa` and two methods implemented here: `mustang_msa`, `mtm_align_msa`.

### opsinmaphmm

```
usage: opsinmaphmm [-h] -i INPUT [-d DATADIR [DATADIR ...]] -o OUTPUT [--pad-n PAD_N] [--pad-c PAD_C] [--max-gap MAX_GAP] [--min-score MIN_SCORE] [-f]
                   [-t THREADS]

Opsin alignments to reference HMM profiles.

options:
  -h, --help            show this help message and exit

Arguments:
  -i INPUT              query sequences in fasta format
  -d DATADIR [DATADIR ...]
                        opsin data directories
  -o OUTPUT             output directory
  --pad-n PAD_N         N-terminal padding for query trimming (default: 30)
  --pad-c PAD_C         C-terminal padding for query trimming (default: 30)
  --max-gap MAX_GAP     Maximum gap for chaining (default: 200)
  --min-score MIN_SCORE
                        Minimum score (default: 15)
  -f                    whether to overwrite files in the output directory if it exists
  -t THREADS            number of threads to use (default: 1)
```

Opsin homology based on protein profiles. The input is a fasta file with multiple opsin sequences and a set of reference dataset from [opsintools-build](https://github.com/BejaLab/opsintools-build/releases). The output is the mapping between positions in the queries and the positions in the reference. For each query the best-fitting dataset is chosen. Both the CLI and the API generate output files located in the specified directory:

* `{dataset}/hmmsearch.txt` - raw `hmmsearch` output for each reference dataset
* `trimmed.pdb` - trimmed query structure
* `opsinmap.json` - json file with the position mapping

The API function `opsinmaphmm` returns the dictionary mapping reference positions to the query positions as its output.

See `opsinmaphmm -h` or `from opsintools import opsinmaphmm; help(opsinmaphmm)` for more details.

### opsinpdb

```
usage: opsinpdb [-h] [-i FILE] [-a FILE] [-o OUTPUT] [-c CHAINS] [-L] [-H] [-W] [-A] [--ligands LIGANDS] [--ligands-remove LIGANDS_REMOVE]

Read a structure or fetch from PDB, cleanup and save as a PDB file.

options:
  -h, --help            show this help message and exit

Arguments:
  -i FILE               Input file
  -a FILE               PDB accession (optionally including chain)
  -o OUTPUT             output path (default: stdout)
  -c CHAINS             comma-separated list of chains (default: all chains)
  -L                    do not remap LYR atoms (default: do remap)
  -H                    do not remove hydrogens (default: remove)
  -W                    do not remove water molecules (default: remove)
  -A                    do not remove alternative conformations (default: remove)
  --ligands LIGANDS     comma-separated list of non-covalent ligands to retain (default: remove all)
  --ligands-remove LIGANDS_REMOVE
                        comma-separated list of covalent ligands to remove (default: don't remove)
```

Fetches a PDB record or reads a local CIF or PDB file, cleans it up and saves as a PDB file. The cleanup includes renaming the atoms of `LYR` (= `LYS+RET`), removing non-covalent ligands, water molecules, hydrogens, alternative states and unwanted chains. See `opsinpdb -h`.
