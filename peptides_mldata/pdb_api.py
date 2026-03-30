"""
RCSB PDB fetching utilities.

Provides functions to:
- fetch_chain_metadata: query RCSB GraphQL API for per-chain entity metadata
  (taxon_id, uniprot_ids) for a given PDB ID.
- fetch_pdb_structure: download the mmCIF file from RCSB and parse all
  amino acid chains (sequence + 3-atom backbone coords).
- fetch_pdb_entry: combine both into a single ready-to-store dict.
"""

import io
import json

import numpy as np
import requests
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser

RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_MMCIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

_CHAIN_METADATA_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    polymer_entities {
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
      }
      rcsb_polymer_entity_container_identifiers {
        uniprot_ids
        auth_asym_ids
      }
    }
  }
}
"""


def fetch_chain_metadata(pdb_id: str) -> dict[str, dict]:
    """
    Fetches per-chain metadata for a single PDB entry.
    Convenience wrapper around fetch_batch_chain_metadata.
    """
    return fetch_batch_chain_metadata([pdb_id]).get(pdb_id.lower(), {})


def fetch_batch_chain_metadata(pdb_ids: list[str], batch_size: int = 200) -> dict[str, dict[str, dict]]:
    """
    Fetches per-chain metadata for a list of PDB IDs in a single GraphQL request.
    Returns {pdb_id: {chain_id: {taxon_id, uniprot_ids}}}.
    Requests are split into batches of `batch_size` to stay within API limits.
    """
    result = {}
    ids_upper = [pid.upper() for pid in pdb_ids]

    for i in range(0, len(ids_upper), batch_size):
        batch = ids_upper[i:i + batch_size]
        response = requests.post(
            RCSB_GRAPHQL_URL,
            json={"query": _CHAIN_METADATA_QUERY, "variables": {"ids": batch}},
            timeout=60,
        )
        response.raise_for_status()
        data = response.json()

        for entry in (data.get("data") or {}).get("entries") or []:
            pdb_id = (entry.get("rcsb_id") or "").lower()
            chain_meta = {}
            for entity in (entry.get("polymer_entities") or []):
                # Collect all unique biological source taxon IDs
                taxon_ids = set()
                source_orgs = entity.get("rcsb_entity_source_organism") or []
                for org in source_orgs:
                    if org.get("ncbi_taxonomy_id"):
                        taxon_ids.add(str(org["ncbi_taxonomy_id"]))

                ids_block = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                uniprot_ids = ids_block.get("uniprot_ids") or []
                uniprot_list = sorted(list(set(i for i in uniprot_ids if i)))

                auth_chains = ids_block.get("auth_asym_ids") or []
                taxon_list = sorted(list(taxon_ids))

                for chain in auth_chains:
                    chain_meta[chain] = {
                        "taxon_ids": taxon_list,
                        "uniprot_ids": uniprot_list,
                    }
            result[pdb_id] = chain_meta

    return result


def fetch_pdb_structure(pdb_id: str) -> dict[str, tuple[str, np.ndarray]]:
    """
    Downloads the mmCIF file from RCSB and parses all amino acid chains.
    Returns {chain_id: (sequence, coords_array)} for standard amino acid chains only.
    Non-protein chains (zero standard residues) are discarded.
    coords_array shape: (N, 3, 3) — N residues × 3 backbone atoms (N, CA, C) × xyz.
    Missing backbone atoms filled with np.nan.
    """
    url = RCSB_MMCIF_URL.format(pdb_id=pdb_id.lower())
    response = requests.get(url, timeout=60)
    response.raise_for_status()

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, io.StringIO(response.text))

    chains = {}
    model = structure[0]

    for chain in model:
        sequence = []
        coords = []

        for residue in chain:
            if residue.id[0] != " ":  # skip HETATM and water
                continue
            resname = residue.resname.capitalize()
            aa = protein_letters_3to1.get(resname)
            if aa is None:
                continue  # skip non-standard residues

            n = residue["N"].coord if "N" in residue else np.full(3, np.nan)
            ca = residue["CA"].coord if "CA" in residue else np.full(3, np.nan)
            c = residue["C"].coord if "C" in residue else np.full(3, np.nan)

            sequence.append(aa)
            coords.append([n, ca, c])

        if not sequence:
            continue  # discard non-amino-acid chains

        seq_str = "".join(sequence)
        coords_arr = np.array(coords, dtype=np.float32)  # shape (N, 3, 3)
        chains[chain.id] = (seq_str, coords_arr)

    return chains


def fetch_pdb_entry(pdb_id: str) -> dict:
    """
    Combines fetch_chain_metadata and fetch_pdb_structure.
    Returns the full entry dict ready for LMDB storage.
    Chains present in structure but absent from metadata get taxon_id=None, uniprot_ids=[].
    """
    pdb_id = pdb_id.lower()
    meta = fetch_chain_metadata(pdb_id)
    structure = fetch_pdb_structure(pdb_id)

    chains = {}
    for chain_id, (sequence, coords) in structure.items():
        m = meta.get(chain_id, {})
        chains[chain_id] = {
            "taxon_id": m.get("taxon_id"),
            "uniprot_ids": m.get("uniprot_ids", []),
            "sequence": sequence,
            "length": len(sequence),
            # Store as list of lists for JSON serialisation; consumer should reconstruct np.array
            "3d_coordinates": coords.tolist(),
        }

    return {
        "pdb_id": pdb_id,
        "chains": chains,
    }
