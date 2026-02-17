def process_pdb(
    file_or_accession: str,
    is_file: bool,
    output_file: str | None,
    chains: list[str],
    atom_map: dict | None,
    non_cov_lig: list[str],
    cov_lig_to_remove: list[str],
    dont_remove_w: bool,
    dont_remove_h: bool,
    dont_remove_alt: bool,
    no_remap: bool):

    import json, re
    from tempfile import TemporaryDirectory
    from pathlib import Path
    import os
    from opsintools.classes.PyMol import PyMol

    def check_alphanum(string):
        return re.match(r'^\w+$', string)

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

    lig_atoms_to_remove = [ "hetatm and not (byres bound_to polymer)" ]
    for lig in non_cov_lig:
        if not check_alphanum(lig):
            raise ValueError(f"Ligand names must be composed only of alphanumerical characters, got: {lig}")
        lig_atoms_to_remove.append(f"not (resn {lig})")

    cov_lig_atoms_to_remove = []
    for lig in cov_lig_to_remove:
        if not check_alphanum(lig):
            raise ValueError(f"Ligand names must be composed only of alphanumerical characters, got: {lig}")
        cov_lig_atoms_to_remove.append(f"(resn {lig} and (byres bound_to polymer))")

    my_space = {}
    with PyMol() as pm:
        if is_file:
            pm.load(file_or_accession, "protein")
        else:
            pm.fetch(file_or_accession, "protein")

        # Remove water and hydrogen atoms
        if not dont_remove_w:
            pm.remove("(solvent)")
        if not dont_remove_h:
            pm.remove("(hydrogens)")
        # Remove alternative states
        if not dont_remove_alt:
            pm.remove("not alt ''+A")
            pm.alter("(all)", "segi, alt = '', ''")
        # remove residues before the start of the protein
        pm.remove("not resi 1-")

        if remove_chains:
            pm.remove(remove_chains)

        # Rename the atoms of the ligand (and of the residue connected to it)
        for resn_from, components in atom_map.items():
            for resn_to, res_data in components.items():
                res_type = res_data['type']
                for atom_from, atom_to in res_data['atoms'].items():
                    my_space['resn_name_type'] = resn_to, atom_to, res_type
                    pm.alter(f"resn {resn_from} and name {atom_from}", "resn, name, type = resn_name_type", space = my_space)
            my_space['extra_atoms'] = []
            pm.iterate(f"resn {resn_from}", "extra_atoms.append(name)", space = my_space)
            if my_space['extra_atoms']:
                raise ValueError(f"Found additional atoms for ligand {resn_from}: {my_space['extra_atoms']}")

        # Remove all extra non-covalent ligands
        pm.remove(' AND '.join(lig_atoms_to_remove))

        # Remove requested covalent ligands
        pm.remove(' OR '.join(cov_lig_atoms_to_remove))

        # Re-number non-polymer resi's
        my_space['polymer_residues'] = set()
        my_space['non_polymer_residues'] = set()

        pm.iterate("polymer", "polymer_residues.add(int(resi))", space = my_space)
        pm.iterate("not polymer", "non_polymer_residues.add(int(resi))", space = my_space)

        if not my_space['polymer_residues']:
            raise ValueError(f"No atoms left, check your filtering criteria")

        my_space['target_resi'] = max(my_space['polymer_residues'])

        # temporally re-number them to put outside of the target range
        my_space['offset'] = my_space['target_resi'] + len(my_space['non_polymer_residues']) + 1000
        pm.alter("not polymer", "resi = offset + int(resi)", space = my_space)
        for resi in my_space['non_polymer_residues']:
            from_resi = resi + my_space['offset']
            my_space['target_resi'] += 1
            pm.alter(f"(not polymer) and resi {from_resi}", "resi = target_resi", space = my_space)

        pm.save(output_file)
