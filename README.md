# opsintools

## opsinmap3d

Opsin homology based on protein structures. The input is an opsin structure in PDB format, the output is the mapping between positions in the query and the positions in the reference in the transmembrane regions. The output directory will contain:

* `aln_to_ref.txt` - alignment of the query to the reference
* `trimmed.pdb` - trimmed query structure
* `trimmed.fasta` - trimmed query sequence
* `t_coffee.fasta` - structural alignment
* `t_coffee.log` - structural alignment log
* `opsinmap.json` - json file with the position mapping

There is an API and a CLI (check out `opsinmap3d --help`).
