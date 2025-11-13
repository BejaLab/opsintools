from tempfile import TemporaryDirectory
import os, re
from pathlib import Path

def check_alphanum(string):
    return re.match(r'^\w+$', string)

def process_pdb(
    file_or_accession: str,
    is_file: bool,
    output_file: str | None,
    chains: list[str],
    atom_map_file: str | None,
    non_cov_lig: list[str],
    dont_remove_w: bool,
    dont_remove_h: bool,
    dont_remove_alt: bool,
    no_remap: bool):

    import json, re
    from psico import exporting
    from pymol import cmd

    cwd = os.getcwd()
    if is_file:
        file_or_accession = Path(file_or_accession).resolve()
    if output_file:
        output_file = Path(output_file).resolve()

    remove_chains = ''
    if chains:
        select = []
        for chain in chains:
            if not check_alphanum(chain):
                raise ValueError(f"Chains must be alphanumeric, got: {chain}")
            select.append(f"not chain {chain}")
        remove_chains = ' and '.join(select)

    remove_hetatms = "hetatm and not (byres bound_to polymer)"
    for lig in non_cov_lig:
        if not check_alphanum(lig):
            raise ValueError(f"Ligand names must be composed only of alphanumerical characters, got: {lig}")
        remove_hetatms += f" and not (resn {lig})"

    atom_map = {}
    if atom_map_file:
        with open(atom_map_file) as file:
            atom_map = json.load(file)

    my_space = {}
    with TemporaryDirectory() as work_dir:
        os.chdir(work_dir)
        if is_file:
            cmd.load(file_or_accession)
        else:
            cmd.fetch(file_or_accession)

        # Remove water and hydrogen atoms
        if not dont_remove_w:
            cmd.remove("(solvent)")
        if not dont_remove_h:
            cmd.remove("(hydrogens)")
        # Remove alternative states
        if not dont_remove_alt:
            cmd.remove("not alt ''+A")
            cmd.alter("(all)", "segi, alt = '', ''")
        # remove residues before the start of the protein
        cmd.remove("not resi 1-")

        if remove_chains:
            cmd.remove(remove_chains)

        # Rename the atoms of the ligand (and of the residue connected to it)
        for resn_from, components in atom_map.items():
            for resn_to, res_data in components.items():
                res_type = res_data['type']
                for atom_from, atom_to in res_data['atoms'].items():
                    my_space['resn_name_type'] = resn_to, atom_to, res_type
                    cmd.alter(f"resn {resn_from} and name {atom_from}", "resn, name, type = resn_name_type", space = my_space)
            my_space['extra_atoms'] = []
            cmd.iterate(f"resn {resn_from}", "extra_atoms.append(name)", space = my_space)
            if my_space['extra_atoms']:
                raise ValueError(f"Found additional atoms for ligand {resn_from}: {my_space['extra_atoms']}")

        # Remove all extra non-covalent ligands
        cmd.remove(remove_hetatms)

        # Re-number heteroatom resi's
        my_space['atom_residues'] = set()
        my_space['hetatm_residues'] = set()

        cmd.iterate("not hetatm", "atom_residues.add(int(resi))", space = my_space)
        cmd.iterate("hetatm", "hetatm_residues.add(int(resi))", space = my_space)

        if not my_space['atom_residues']:
            raise ValueError(f"No atoms left, check your filtering criteria")

        my_space['target_resi'] = max(my_space['atom_residues'])

        # temporally re-number them to put outside of the target range
        my_space['offset'] = my_space['target_resi'] + len(my_space['hetatm_residues']) + 1000
        cmd.alter("hetatm", "resi = offset + int(resi)", space = my_space)
        for resi in my_space['hetatm_residues']:
            from_resi = resi + my_space['offset']
            my_space['target_resi'] += 1
            cmd.alter(f"hetatm and resi {from_resi}", "resi = target_resi", space = my_space)

        cmd.sort()
        cmd.set('pdb_conect_nodup', 0)
        if output_file:
            exporting.save_pdb(output_file, seqres = True)
        else:
            exporting.save_pdb('output.pdb', seqres = True)
            with open('output.pdb') as file:
                print(file.read(), end = '')
