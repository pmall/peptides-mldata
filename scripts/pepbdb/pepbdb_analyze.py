import argparse
import numpy as np
from peptides_mldata.iterators.pepbdb import iter_pepbdb

# High-quality filtering parameters
RESOL_MIN = 0.001
RESOL_MAX = 2.5
MOL_TYPE = "prot"
NONSTANDARD = False

SITE_BIN_SIZE = 32

def analyze_pepbdb(min_pep_len: int, max_pep_len: int):
    print(f"High-Quality PepBDB Analysis")
    print(f"Filters: Pep [{min_pep_len}-{max_pep_len}] | Res [{RESOL_MIN}-{RESOL_MAX}] | Type: {MOL_TYPE} | Std AA: {not NONSTANDARD}\n")

    site_lengths = []
    total_yielded = 0

    # Using the standard iterator filters
    iterator = iter_pepbdb(
        archive_path="data/pepbdb-20200318.zip",
        nonstandard_aa=NONSTANDARD,
        resolution_min=RESOL_MIN,
        resolution_max=RESOL_MAX,
        mol_type=MOL_TYPE,
        min_pep_len=min_pep_len,
        max_pep_len=max_pep_len,
        verbose=False # Keep it clean since we are analyzing
    )

    for item in iterator:
        total_yielded += 1
        site_lengths.append(item['target']['length'])

        if total_yielded % 500 == 0:
            print(f"  {total_yielded} matches found...")

    print(f"\nTotal high-quality matches found: {len(site_lengths)}\n")

    if not site_lengths:
        print("No matches found within these strict quality constraints.")
        return

    data = np.array(site_lengths)
    max_site = int(np.max(data))
    bins = list(range(0, max_site + SITE_BIN_SIZE + 1, SITE_BIN_SIZE))
    counts, edges = np.histogram(data, bins=bins)
    max_count = max(counts) if max(counts) > 0 else 1
    
    print(f"{'Binder Bin':<12} | {'Count':>5} | {'%':>4} | Visualization")
    print("-" * 65)

    for i in range(len(counts)):
        # Skip empty tail logic
        if counts[i] == 0 and i > 12: 
            if not any(counts[i:]):
                break
        
        blo = int(edges[i])
        bhi = int(edges[i + 1])
        percent = (counts[i] / len(site_lengths)) * 100
        hashes = "#" * int(counts[i] / max_count * 40)
        
        print(f"[{blo:>4}-{bhi:<4}) | {counts[i]:>5} | {percent:>3.0f}% | {hashes}")

    print("\nTotal High-Quality Sample Count:", len(site_lengths))
    print(f"Mean: {np.mean(data):.1f} | Median: {np.median(data):.0f} | Max: {max_site}")


def main():
    parser = argparse.ArgumentParser(description="Analyze PepBDB binding site lengths within peptide length constraints.")
    parser.add_argument("--min_pep_len", type=int, default=4, help="Minimum peptide length (default: 4)")
    parser.add_argument("--max_pep_len", type=int, default=32, help="Maximum peptide length (default: 32)")
    args = parser.parse_args()

    analyze_pepbdb(args.min_pep_len, args.max_pep_len)

if __name__ == "__main__":
    main()
