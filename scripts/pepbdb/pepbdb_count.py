import argparse
from peptides_mldata.iterators.pepbdb import iter_pepbdb

RESOLUTION_MIN = 0.001
RESOLUTION_MAX = 2.5
MOL_TYPE = "prot"
NONSTANDARD = False


def count_pepbdb(min_pep_len: int, max_pep_len: int, min_target_len: int, max_target_len: int, archive_path: str):
    print(f"PepBDB Sweet Spot Count: Pep [{min_pep_len}-{max_pep_len}] | Site [{min_target_len}-{max_target_len}]\n")

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
    )

    total = 0
    passed = 0
    for item in iterator:
        total += 1
        passed += 1
        if total % 1000 == 0:
            print(f"  {total} entries scanned, {passed} matches...")

    print(f"\nTotal scanned: {total}")
    print(f"Sweet spot entries: {passed} ({passed/total*100:.1f}%)" if total > 0 else "\nNo matches found.")


def main():
    parser = argparse.ArgumentParser(description="Count PepBDB entries within length constraints.")
    parser.add_argument("archive_path", nargs="?", default="data/pepbdb-20200318.zip",
                        help="Path to the PepBDB zip archive (default: data/pepbdb-20200318.zip)")
    parser.add_argument("--min_pep_len", type=int, default=4, help="Minimum peptide length (default: 4)")
    parser.add_argument("--max_pep_len", type=int, default=32, help="Maximum peptide length (default: 32)")
    parser.add_argument("--min_target_len", type=int, default=15, help="Minimum target/site length (default: 15)")
    parser.add_argument("--max_target_len", type=int, default=280, help="Maximum target/site length (default: 280)")
    args = parser.parse_args()

    count_pepbdb(args.min_pep_len, args.max_pep_len, args.min_target_len, args.max_target_len, args.archive_path)


if __name__ == "__main__":
    main()
