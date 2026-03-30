"""
Batch download full PDB entries from RCSB and store them in an LMDB database.

Input : any CSV file with pdb_id as the first column (header expected, duplicates ok).
Output: LMDB database where each key is a pdb_id (bytes) and each value is
        a JSON-encoded dict with all amino acid chains, their sequences,
        3-atom backbone coordinates, and entity metadata (taxon_id, uniprot_ids).

Resume-safe: already-stored keys are skipped on re-run.
"""

import argparse
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import lmdb

from peptides_mldata.pdb_api import fetch_batch_chain_metadata, fetch_pdb_structure

LMDB_MAP_SIZE = 50 * 1024 ** 3  # 50 GB max map (sparse, actual disk usage is minimal)


def fetch_pdb(lmdb_path: str, csv_path: str, workers: int = 8, verbose: bool = False):
    """
    Read unique PDB IDs from the first column of csv_path,
    download and parse each entry via RCSB, and store in lmdb_path.

    Metadata (taxon_id, uniprot_ids) is fetched in a single batched GraphQL request.
    mmCIF structural files are downloaded in parallel with a thread pool.
    """
    # --- Collect unique PDB IDs ---
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header
        pdb_ids = list({row[0].strip().lower() for row in reader if row})

    print(f"Found {len(pdb_ids)} unique PDB IDs in {csv_path}")

    # --- Open LMDB and filter already-stored keys ---
    env = lmdb.open(lmdb_path, map_size=LMDB_MAP_SIZE, max_dbs=0)
    with env.begin() as txn:
        todo = [pid for pid in pdb_ids if txn.get(pid.encode()) is None]
    print(f"Skipping {len(pdb_ids) - len(todo)} already stored. Fetching {len(todo)}...")

    if not todo:
        print("Nothing to do.")
        env.close()
        return

    # --- Batch-fetch all metadata upfront (1 GraphQL request per 200 entries) ---
    print(f"Fetching chain metadata for {len(todo)} entries...")
    all_meta = fetch_batch_chain_metadata(todo)
    print("Metadata ready. Downloading mmCIF structures...")

    # --- Parallel mmCIF download + parse ---
    done = 0
    errors = 0

    def _fetch_structure(pdb_id: str):
        structure = fetch_pdb_structure(pdb_id)
        meta = all_meta.get(pdb_id, {})
        chains = {}
        for chain_id, (sequence, coords) in structure.items():
            m = meta.get(chain_id, {})
            chains[chain_id] = {
                "taxon_ids": m.get("taxon_ids", []),
                "uniprot_ids": m.get("uniprot_ids", []),
                "sequence": sequence,
                "length": len(sequence),
                "3d_coordinates": coords.tolist(),
            }
        return {"pdb_id": pdb_id, "chains": chains}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_fetch_structure, pid): pid for pid in todo}
        for future in as_completed(futures):
            pdb_id = futures[future]
            try:
                entry = future.result()
                value = json.dumps(entry).encode("utf-8")
                with env.begin(write=True) as txn:
                    txn.put(pdb_id.encode(), value)
                done += 1
                if verbose:
                    chains = list(entry["chains"].keys())
                    print(f"  [OK] {pdb_id}: {len(chains)} chains ({', '.join(chains)})")
                elif done % 100 == 0:
                    print(f"  {done}/{len(todo)} done, {errors} errors...")
            except Exception as e:
                errors += 1
                print(f"  [FAIL] {pdb_id}: {e}", file=sys.stderr)

    env.close()
    print(f"\nDone. {done} stored, {errors} failed.")


def main():
    parser = argparse.ArgumentParser(
        description="Batch download PDB entries from RCSB and store in LMDB."
    )
    parser.add_argument("input_csv",
                        help="Path to a CSV file with pdb_id as the first column")
    parser.add_argument("lmdb_path", nargs="?", default="data/pdb.lmdb",
                        help="Path to the output LMDB (default: data/pdb.lmdb)")
    parser.add_argument("--workers", type=int, default=8,
                        help="Number of parallel download threads (default: 8)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print per-entry status")
    args = parser.parse_args()

    fetch_pdb(args.lmdb_path, args.input_csv, args.workers, args.verbose)


if __name__ == "__main__":
    main()
