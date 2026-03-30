"""
Inspect pdb.lmdb and count the distribution of (num_taxons, num_proteins) per chain.
"""

import json
from collections import Counter
import lmdb
import argparse

def inspect_pdb(lmdb_path: str, show_examples: bool = False):
    print(f"Inspecting {lmdb_path}...")
    
    try:
        env = lmdb.open(lmdb_path, readonly=True, lock=False)
    except Exception as e:
        print(f"Error: Could not open {lmdb_path}: {e}")
        return
    
    # Separate counters for taxon_ids and uniprot_ids
    taxon_counts = Counter()
    protein_counts = Counter()
    total_entries = 0
    total_chains = 0
    
    no_taxon_examples = []
    no_protein_examples = []

    with env.begin() as txn:
        cursor = txn.cursor()
        for key, value in cursor:
            total_entries += 1
            try:
                entry = json.loads(value.decode("utf-8"))
                for chain_id, cdata in entry.get("chains", {}).items():
                    total_chains += 1
                    n_tax = len(cdata.get("taxon_ids", []))
                    n_prot = len(cdata.get("uniprot_ids", []))
                    
                    taxon_counts[n_tax] += 1
                    protein_counts[n_prot] += 1
                    
                    if n_tax == 0 and len(no_taxon_examples) < 10:
                        no_taxon_examples.append(f"{entry.get('pdb_id', key.decode())}_{chain_id}")
                    if n_prot == 0 and len(no_protein_examples) < 10:
                        no_protein_examples.append(f"{entry.get('pdb_id', key.decode())}_{chain_id}")
            except Exception:
                continue

    env.close()

    # Display results
    print("\n" + "="*40)
    print(f"Summary: {total_entries:,} entries, {total_chains:,} chains")
    print("="*40)

    print(f"\nTaxon IDs per chain:")
    print("-" * 40)
    print(f"{'Count':<10} | {'Chains':<10} | {'Percentage'}")
    print("-" * 40)
    for val in sorted(taxon_counts.keys()):
        count = taxon_counts[val]
        pct = (count / total_chains) * 100 if total_chains > 0 else 0
        print(f"{val:<10} | {count:10,} | {pct:8.2f}%")

    print(f"\nProtein IDs per chain:")
    print("-" * 40)
    print(f"{'Count':<10} | {'Chains':<10} | {'Percentage'}")
    print("-" * 40)
    for val in sorted(protein_counts.keys()):
        count = protein_counts[val]
        pct = (count / total_chains) * 100 if total_chains > 0 else 0
        print(f"{val:<10} | {count:10,} | {pct:8.2f}%")

    if show_examples:
        print(f"\nExamples of no Taxon IDs: {', '.join(no_taxon_examples)}")
        print(f"Examples of no Protein IDs: {', '.join(no_protein_examples)}")
    print("\n")

def main():
    parser = argparse.ArgumentParser(description="Inspect PDB LMDB metadata distribution.")
    parser.add_argument("lmdb_path", nargs="?", default="data/pdb.lmdb",
                        help="Path to the PDB LMDB (default: data/pdb.lmdb)")
    parser.add_argument("--examples", action="store_true", help="Show example chain IDs with missing metadata")
    args = parser.parse_args()
    inspect_pdb(args.lmdb_path, args.examples)

if __name__ == "__main__":
    main()
