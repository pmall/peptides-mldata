"""
PepBDB data iterator.

Column explanation (verified from peptidelist.txt):
Col 1 [0] - PDB ID
Col 2 [1] - peptide chain ID
Col 3 [2] - peptide length
Col 4 [3] - number of atoms in peptide
Col 5 [4] - protein chain ID
Col 6 [5] - number of atoms in protein
Col 7 [6] - number of atom contacts
Col 8 [7] - peptide with nonstandard amino acid? (1 means yes)
Col 9 [8] - (unused/skipped)
Col 10[9] - resolution (-1.00 means NMR structure)
Col 11[10]- molecular type

Typical source file: data/pepbdb-20200318.zip
"""

import io
import zipfile
import numpy as np
from typing import Iterator, Dict, Any, Tuple
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1


def _parse_pdb_stream(stream, target_chain_id: str) -> Tuple[str, np.ndarray]:
    """
    Parses an open stream of PDB data using BioPython.
    Returns (sequence_str, coords_array).
    """
    parser = PDBParser(QUIET=True)
    content = stream.read().decode("utf-8")
    structure = parser.get_structure("x", io.StringIO(content))
    model = structure[0]
    
    if target_chain_id not in model:
        raise KeyError(f"Chain {target_chain_id} not found in PDB")
        
    chain = model[target_chain_id]
    sequence = []
    coords = []

    for residue in chain:
        if residue.id[0] != ' ':  # skip hetero atoms
            continue

        resname = residue.resname.capitalize()
        aa_char = protein_letters_3to1.get(resname, 'X')

        n_coord = residue['N'].coord if 'N' in residue else np.full(3, np.nan)
        ca_coord = residue['CA'].coord if 'CA' in residue else np.full(3, np.nan)
        c_coord = residue['C'].coord if 'C' in residue else np.full(3, np.nan)

        sequence.append(aa_char)
        coords.append([n_coord, ca_coord, c_coord])

    seq_str = "".join(sequence)
    coords_arr = np.array(coords, dtype=np.float32)
    return seq_str, coords_arr


def iter_pepbdb(
    archive_path: str = "data/pepbdb-20200318.zip",
    nonstandard_aa: bool = None,
    resolution_min: float = None,
    resolution_max: float = None,
    mol_type: str = None,
    min_pep_len: int = None,
    max_pep_len: int = None,
    min_target_len: int = None,
    max_target_len: int = None,
    verbose: bool = False
) -> Iterator[Dict[str, Any]]:
    """
    Yields parsed records from the PepBDB zip archive on the fly.

    Uses zipfile random access: iterates peptidelist.txt line by line,
    applies metadata filters, then fetches receptor.pdb and peptide.pdb
    by name for each passing entry and yields immediately.
    Zero buffering of PDB data.

    Malformed or missing entries are skipped to keep the iterator alive.
    Set verbose=True to print exceptions to stdout.
    """
    with zipfile.ZipFile(archive_path, "r") as zf:
        with zf.open("peptidelist.txt") as f:
            for line in io.TextIOWrapper(f, encoding="utf-8"):
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 11:
                    if verbose:
                        print(f"[pepbdb] Skip malformed index entry: {line}")
                    continue

                pdb_id = parts[0]
                pep_chain = parts[1]
                prot_chain = parts[4]
                entry_nonstandard_aa = int(parts[7])
                resolution = float(parts[9])
                entry_mol_type = parts[10]
                folder = f"{pdb_id}_{pep_chain}"

                # Apply metadata filters before opening any PDB file
                if nonstandard_aa is not None and bool(entry_nonstandard_aa) != nonstandard_aa:
                    continue
                if resolution_min is not None and resolution <= resolution_min:
                    continue
                if resolution_max is not None and resolution > resolution_max:
                    continue
                if mol_type is not None and entry_mol_type != mol_type:
                    continue

                if min_pep_len is not None and int(parts[2]) < min_pep_len:
                    continue
                if max_pep_len is not None and int(parts[2]) > max_pep_len:
                    continue

                try:
                    receptor_path = f"pepbdb/{folder}/receptor.pdb"
                    peptide_path = f"pepbdb/{folder}/peptide.pdb"

                    with zf.open(receptor_path) as rec_stream:
                        rec_seq, rec_coords = _parse_pdb_stream(rec_stream, prot_chain)

                    if min_target_len is not None and len(rec_seq) < min_target_len:
                        continue
                    if max_target_len is not None and len(rec_seq) > max_target_len:
                        continue

                    with zf.open(peptide_path) as pep_stream:
                        pep_seq, pep_coords = _parse_pdb_stream(pep_stream, pep_chain)
                except Exception as e:
                    if verbose:
                        print(f"[pepbdb] Skip {folder}: {e}")
                    continue

                yield {
                    "source_key": folder,
                    "source_db": "pepbdb",
                    "pdb_id": pdb_id,
                    "target": {
                        "chain": prot_chain,
                        "length": len(rec_seq),
                        "sequence": rec_seq,
                        "3d_coordinates": rec_coords
                    },
                    "peptide": {
                        "chain": pep_chain,
                        "length": len(pep_seq),
                        "sequence": pep_seq,
                        "3d_coordinates": pep_coords
                    }
                }
