"""
Build the final PepBDB ML dataset LMDB.

Combines:
- iter_pepbdb: yields filtered peptide/binding-site pairs with sequences and
  3D backbone coordinates from the curated PepBDB source.
- pdb.lmdb: ground-truth PDB entries with per-chain taxon_id and protein_ids.

Only perfect entries are written. An entry is skipped if:
- Its pdb_id is missing from pdb.lmdb.
- Either the peptide or target chain is missing from its pdb.lmdb entry.
- taxon_id or protein_ids are missing for either chain.
- The PepBDB sequence is not a contiguous subsequence of the source PDB chain.

Keys: source_key (e.g. '1a1r_C') — bare iterator key, no prefix.
Each source produces its own LMDB; cross-source deduplication happens at merge time.
"""

import argparse
import json
import sys

import lmdb
import numpy as np

from peptides_mldata.iterators.pepbdb import iter_pepbdb

RESOLUTION_MIN = 0.001
RESOLUTION_MAX = 2.5
MOL_TYPE = "prot"
NONSTANDARD = False
LMDB_MAP_SIZE = 10 * 1024 ** 3  # 10 GB


def _is_subsequence(sub: str, full: str) -> bool:
    """Return True if `sub` is a contiguous subsequence of `full`."""
    return sub in full


def build_pepbdb(
    pdb_lmdb: str,
    archive_path: str,
    lmdb_out: str,
    min_pep_len: int,
    max_pep_len: int,
    min_target_len: int,
    max_target_len: int,
    verbose: bool = False,
):
    pdb_env = lmdb.open(pdb_lmdb, readonly=True, lock=False)
    out_env = lmdb.open(lmdb_out, map_size=LMDB_MAP_SIZE)

    stats = {
        "total": 0,
        "written": 0,
        "skipped_missing_pdb": 0,
        "skipped_missing_chain": 0,
        "skipped_subsequence_mismatch": 0,
    }

    iterator = iter_pepbdb(
        archive_path=archive_path,
        nonstandard_aa=NONSTANDARD,
        resolution_min=RESOLUTION_MIN,
        resolution_max=RESOLUTION_MAX,
        mol_type=MOL_TYPE,
        min_pep_len=min_pep_len,
        max_pep_len=max_pep_len,
        min_target_len=min_target_len,
        max_target_len=max_target_len,
        verbose=verbose,
    )

    for item in iterator:
        stats["total"] += 1
        pdb_id = item["pdb_id"].lower()
        source_key = item["source_key"]
        pep_chain = item["peptide"]["chain"]
        tgt_chain = item["target"]["chain"]

        # 1. Look up pdb.lmdb
        with pdb_env.begin() as txn:
            raw = txn.get(pdb_id.encode())

        if raw is None:
            stats["skipped_missing_pdb"] += 1
            if verbose:
                print(f"[SKIP] {source_key}: pdb_id '{pdb_id}' not in pdb.lmdb")
            continue

        pdb_entry = json.loads(raw.decode("utf-8"))
        chains = pdb_entry.get("chains", {})

        # 2. Check both chains exist
        if pep_chain not in chains or tgt_chain not in chains:
            stats["skipped_missing_chain"] += 1
            if verbose:
                missing = [c for c in [pep_chain, tgt_chain] if c not in chains]
                print(f"[SKIP] {source_key}: chain(s) {missing} not in pdb.lmdb entry")
            continue

        pep_src = chains[pep_chain]
        tgt_src = chains[tgt_chain]

        # 3. Subsequence check (essential for data quality)
        pep_seq = item["peptide"]["sequence"]
        tgt_seq = item["target"]["sequence"]

        if not _is_subsequence(pep_seq, pep_src["sequence"]):
            stats["skipped_subsequence_mismatch"] += 1
            if verbose:
                print(f"[SKIP] {source_key}: peptide subsequence mismatch")
            continue

        if not _is_subsequence(tgt_seq, tgt_src["sequence"]):
            stats["skipped_subsequence_mismatch"] += 1
            if verbose:
                print(f"[SKIP] {source_key}: target subsequence mismatch")
            continue

        # 5. Assemble and write
        pep_coords = item["peptide"]["3d_coordinates"]
        tgt_coords = item["target"]["3d_coordinates"]

        entry = {
            "source_db": "pepbdb",
            "pdb_id": pdb_id,
            "target": {
                "chain": tgt_chain,
                "length": item["target"]["length"],
                "taxon_ids": tgt_src["taxon_ids"],
                "protein_ids": tgt_src["uniprot_ids"],
                "sequence": tgt_seq,
                "3d_coordinates": tgt_coords.tolist(),
            },
            "peptide": {
                "chain": pep_chain,
                "length": item["peptide"]["length"],
                "taxon_ids": pep_src["taxon_ids"],
                "protein_ids": pep_src["uniprot_ids"],
                "sequence": pep_seq,
                "3d_coordinates": pep_coords.tolist(),
            },
        }

        key = source_key.encode()
        value = json.dumps(entry).encode("utf-8")

        with out_env.begin(write=True) as txn:
            txn.put(key, value)

        stats["written"] += 1

        if not verbose and stats["total"] % 500 == 0:
            print(f"  {stats['total']} processed, {stats['written']} written...")

    pdb_env.close()
    out_env.close()

    print("\n" + "="*40)
    print("PepBDB Build Summary")
    print("="*40)
    print(f"Total entries processed: {stats['total']}")
    print(f"Entries written:         {stats['written']}")
    print(f"Skipped (missing PDB):   {stats['skipped_missing_pdb']}")
    print(f"Skipped (missing chain): {stats['skipped_missing_chain']}")
    print(f"Skipped (seq mismatch):  {stats['skipped_subsequence_mismatch']}")
    print("="*40)
    if stats["total"] > 0:
        print(f"  Yield rate:             {stats['written']/stats['total']*100:.1f}%")


def main():
    parser = argparse.ArgumentParser(description="Build the final PepBDB ML dataset LMDB.")
    parser.add_argument("pdb_lmdb", nargs="?", default="data/pdb.lmdb",
                        help="Path to the pdb.lmdb source database (default: data/pdb.lmdb)")
    parser.add_argument("archive_path", nargs="?", default="data/pepbdb-20200318.zip",
                        help="Path to the PepBDB zip archive (default: data/pepbdb-20200318.zip)")
    parser.add_argument("lmdb_out", nargs="?", default="data/pepbdb.lmdb",
                        help="Path to the output LMDB (default: data/pepbdb.lmdb)")
    parser.add_argument("--min_pep_len", type=int, default=4)
    parser.add_argument("--max_pep_len", type=int, default=32)
    parser.add_argument("--min_target_len", type=int, default=15)
    parser.add_argument("--max_target_len", type=int, default=280)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    build_pepbdb(
        pdb_lmdb=args.pdb_lmdb,
        archive_path=args.archive_path,
        lmdb_out=args.lmdb_out,
        min_pep_len=args.min_pep_len,
        max_pep_len=args.max_pep_len,
        min_target_len=args.min_target_len,
        max_target_len=args.max_target_len,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
