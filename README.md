# opsintools

This project contains a set of tools for analysis of opsin structures and sequences.

Each tool has an API and a CLI. Keep in mind that the interfaces might change in the future since the project is still under active development.

# Installation

There are python and non-python dependencies: the most straightforward installation is via conda which takes care of both.

## conda

Use the provided environment definition: `conda env create -n opsintools -f env.yaml` and then `conda activate opsintools`.

## pipy

`opsintools` and its python dependencies can be installed with `pip install git+https://github.com/BejaLab/opsintools`. Non-python dependencies should then be installed separately. As of now these include:

* [`dssp`](https://swift.cmbi.umcn.nl/gv/dssp/index.html) (v.3.0.0 tested)
* [`t-coffee`](https://tcoffee.crg.eu/) (v.13.46 tested)
* [`US-align`](https://github.com/pylelab/USalign) (v.20241201 test)

## opsinmap3d

Opsin homology based on protein structures. The input is an opsin structure in PDB format and a suitable reference dataset (can be dowloaded from [opsintools-build](https://github.com/BejaLab/opsintools-build/releases)), the output is the mapping between positions in the query and the positions in the reference in the transmembrane regions. The output directory will contain:

* `aln_to_ref.txt` - alignment of the query to the reference
* `trimmed.pdb` - trimmed query structure
* `t_coffee.aln` - structural alignment
* `t_coffee.aln.score_ascii` - alignment scores
* `t_coffee.aln.log` - structural alignment log
* `opsinmap.json` - json file with the position mapping

## opsinalign3d

A user-friendly wrapper for running t-coffee on a set of PDB files and parsing the results. The output is an object of the custom class Tcoffee. The output directory will contain:

* `t_coffee.aln` - structural alignment
* `t_coffee.log` - structural alignment log
* `opsinmap.json` - json file with the position mapping
