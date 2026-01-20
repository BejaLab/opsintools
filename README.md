# opsintools

This project contains a set of tools for analysis of opsin structures and sequences.

Each tool has an API and a CLI. Keep in mind that the interfaces might change in the future since the project is under active development.

# Installation

There are python and non-python dependencies: the most straightforward installation is via conda which takes care of both.

## conda

Use the provided environment definition: `conda env create -n opsintools -f env.yaml --channel-priority flexible` and then `conda activate opsintools`. Use `mamba` if available.

## pypi

`opsintools` and its python dependencies can be installed with `pip install git+https://github.com/BejaLab/opsintools`. Non-python dependencies should then be installed separately. As of now these include:

* [`dssp`](https://swift.cmbi.umcn.nl/gv/dssp/index.html) (v.3.0.0 tested)
* [`t-coffee`](https://tcoffee.crg.eu/) (v.13.46 tested)
* [`US-align`](https://github.com/pylelab/USalign) (v.20241201 tested)

## opsinmap3d

Opsin homology based on protein structures. The input is an opsin structure in PDB format and a suitable reference dataset (can be dowloaded from [opsintools-build](https://github.com/BejaLab/opsintools-build/releases)), the output is the mapping between positions in the query and the positions in the reference in the transmembrane regions. Both the CLI and the API generate output files located in the specified directory:

* `aln_to_ref.txt` - alignment of the query to the reference
* `trimmed.pdb` - trimmed query structure
* `t_coffee.aln` - structural alignment
* `t_coffee.aln.score_ascii` - alignment scores
* `t_coffee.aln.log` - structural alignment log
* `opsinmap.json` - json file with the position mapping

The API function `opsinmap3d` returns the dictionary mapping reference positions to the query positions as its output.

See `opsinmap3d -h` or `from opsintools import opsinmap3d; help(opsinmap3d)` for more details.

## opsinalign3d

Runs `t-coffee` on a set of PDB files and parses the results. The output is an object of the custom class `Tcoffee`. The output directory will contain:

* `t_coffee.aln` - structural alignment
* `t_coffee.log` - structural alignment log

See `opsinalign3d -h` or `from opsintools import opsinalign3d; help(opsinalign3d)` for more details.

Currently supported methods are: `sap_pair`, `mustang_pair`, `t_coffee_msa`, `probcons_msa`, `vmaffteinsi_msa`, `mafftfftns1_msa`, `mafftfftnsi_msa`, `mafftginsi_msa`, `mafftlinsi_msa`, `mafft_msa`, `mafftnwnsi_msa`, `mafftsparsecore_msa` and two methods implemented here: `mustang_msa`, `mtm_align_msa`.

## opsinmaphmm

Opsin homology based on protein profiles. The input is a fasta file with multiple opsin sequences and a set of reference dataset from [opsintools-build](https://github.com/BejaLab/opsintools-build/releases). The output is the mapping between positions in the queries and the positions in the reference. For each query the best-fitting dataset is chosen. Both the CLI and the API generate output files located in the specified directory:

* `{dataset}/hmmsearch.txt` - raw `hmmsearch` output for each reference dataset
* `trimmed.pdb` - trimmed query structure
* `opsinmap.json` - json file with the position mapping

The API function `opsinmaphmm` returns the dictionary mapping reference positions to the query positions as its output.

See `opsinmaphmm -h` or `from opsintools import opsinmaphmm; help(opsinmaphmm)` for more details.

## opsinpdb

Fetches a PDB record or reads a local CIF or PDB file, cleans it up and saves as a PDB file. The cleanup includes renaming the atoms of `LYR` (= `LYS+RET`), removing non-covalent ligands, water molecules, hydrogens, alternative states and unwanted chains. See `opsinpdb -h`.
