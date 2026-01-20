from .globals import run_with_logger, logger, create_output_dir, read_database
from .globals import BUILD_REPO

def opsindata(output_dir: str, repo: str = BUILD_REPO, tag: str = None, force: bool = False):

    import argparse
    import os
    import requests
    import zipfile
    from pathlib import Path

    # GitHub API
    base_url = f"https://api.github.com/repos/{repo}/releases"
    suffix = 'latest' if tag is None or tag == 'latest' else f"tags/{tag}"
    release_url = f"{base_url}/{suffix}"

    resp = requests.get(release_url, headers={"Accept": "application/vnd.github.v3+json"})
    if resp.status_code != 200:
        raise SystemExit(f"Error fetching release: {resp.status_code} – {resp.reason}\n{resp.text}")

    release = resp.json()
    tag_name = release['tag_name']
    logger.info(f"Target release: {tag_name}")

    assets = release.get('assets', [])
    if not assets:
        logger.error("No assets found")
        return

    subdir = Path(output_dir) / str(tag_name)

    if subdir.is_file():
        logger.error(f"{subdir} already exists and is a file")
        return
    if subdir.is_dir() and any(subdir.iterdir()) and not force:
        logger.info(f"Folder '{subdir}' already exists and contains files - skipping (use --force to override)")
        return
    subdir.mkdir(parents = True, exist_ok = True)

    for asset in assets:
        name = asset['name']
        local_path = subdir / name

        logger.info(f"Downloading {name}")
        r = requests.get(asset['browser_download_url'], stream=True)
        r.raise_for_status()
        with open(local_path, 'wb') as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        logger.info(f"Downloaded {name}")

        # Unpack & clean only if it's a zip
        if name.lower().endswith('.zip'):
            try:
                with zipfile.ZipFile(local_path, 'r') as z:
                    z.extractall(subdir)
            except zipfile.BadZipFile:
                logger.error(f"{name} is not a valid zip")
                return
            except Exception as e:
                logger.error(f"Failed to unpack/remove {name}: {e}")
                return
            finally:
                local_path.unlink()

    logger.info(f"Finished processing release {tag_name}")

def opsindata_cli():
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Fetch the reference datasets.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-d', metavar = 'DIRNAME', required = True, help = 'download directory')
    main_group.add_argument('-r', metavar = 'REPO', type = str, default = BUILD_REPO, help = 'custom gh repository')
    main_group.add_argument('-t', metavar = 'REPO', type = str, default = 'latest', help = 'release tag name')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to force download')

    args = parser.parse_args()

    run_with_logger(opsindata,
        output_dir = args.d,
        tag = args.t,
        repo = args.r,
        force = args.f
    )

def opsinpdb(
    file_or_accession: str,
    output_file: str,
    is_file: bool,
    chains: list[str] = [],
    non_cov_lig: list[str] = [],
    dont_remove_w: bool = False,
    dont_remove_h: bool = False,
    dont_remove_alt: bool = False,
    no_remap: bool = False) -> None:
    """Read a structure file or fetch from PDB, cleaunp and save as PDB file

    :param accession: PDB accession
    :param output_file: output file path
    :no_remap: whether to remap LYR atoms
    """
    from pathlib import Path

    if not file_or_accession:
        raise ValueError("No file or accession provided")
    if is_file and not Path(file_or_accession).is_file():
        raise ValueError(f"File {file_or_accession} does not exist")
    
    import importlib.resources
    from opsintools.scripts.pdb import process_pdb

    atom_map_file = None if no_remap else importlib.resources.files('opsintools.resources').joinpath('atom_map.json')

    process_pdb(
        file_or_accession = file_or_accession,
        is_file = is_file,
        output_file = output_file,
        chains = chains,
        atom_map_file = atom_map_file,
        non_cov_lig = non_cov_lig,
        dont_remove_w = dont_remove_w,
        dont_remove_h = dont_remove_h,
        dont_remove_alt = dont_remove_alt,
        no_remap = no_remap)

def opsinpdb_cli():
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Read a structure or fetch from PDB, cleanup and save as a PDB file.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', metavar = 'FILE', help = 'Input file')
    main_group.add_argument('-a', metavar = 'FILE', help = 'PDB accession (optionally including chain)')

    main_group.add_argument('-o', metavar = 'OUTPUT', required = False, help = 'output path (default: stdout)')
    main_group.add_argument('-c', metavar = 'CHAINS', required = False, help = 'comma-separated list of chains (default: all chains)')
    main_group.add_argument('-L', action = 'store_true', help = 'do not remap LYR atoms (default: do remap)')
    main_group.add_argument('-H', action = 'store_true', help = 'do not remove hydrogens (default: remove)')
    main_group.add_argument('-W', action = 'store_true', help = 'do not remove water molecules (default: remove)')
    main_group.add_argument('-A', action = 'store_true', help = 'do not remove alternative conformations (default: remove)')
    main_group.add_argument('--ligands', metavar = 'LIGANDS', help = 'comma-separated list of non-covalent ligands to retain (default: remove all)')

    args = parser.parse_args()

    chains = [] if not args.c else args.c.split(',')
    ligands = [] if not args.ligands else args.ligands.split(',')
    if bool(args.i) == bool(args.a):
        raise ValueError("Must specify either -a or -i")

    if args.i:
        is_file = True
        file_or_accession = args.i
    else:
        is_file = False
        file_or_accession = args.a
    output_file = args.o if args.o else '/dev/stdout'

    run_with_logger(opsinpdb,
        file_or_accession = file_or_accession,
        output_file = args.o,
        is_file = is_file,
        chains = chains,
        no_remap = args.L,
        non_cov_lig = ligands,
        dont_remove_w = args.W,
        dont_remove_h = args.H,
        dont_remove_alt = args.A
    )
