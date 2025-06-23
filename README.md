# opsintools

## opsinmap3d

Opsin homology based on protein structures. The input is an opsin structure in PDB format, the output is a json file (`opsinmap.json`) mapping the positions in the query to the positions in the reference in the transmembrane regions. Additional output files are:

* `aln_to_ref.txt` - alignment of the query to the reference
* `trimmed.pdb` - trimmed query structure
* `trimmed.fasta` - trimmed query sequence
* `t_coffee.fasta` - structural alignment
* `t_coffee.log` - structural alignment log
