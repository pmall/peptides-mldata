import zipfile
import io
import argparse
from peptides_mldata.iterators.pepbdb import iter_pepbdb


def count_peptidelist_lines(zip_path):
    count = 0
    with zipfile.ZipFile(zip_path, 'r') as zf:
        with zf.open("peptidelist.txt") as f:
            for line in io.TextIOWrapper(f, encoding="utf-8"):
                if line.strip():
                    count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description="Verify PepBDB iterator yields all entries from the index.")
    parser.add_argument("archive_path", nargs="?", default="data/pepbdb-20200318.zip",
                        help="Path to the PepBDB zip archive (default: data/pepbdb-20200318.zip)")
    parser.add_argument("--verbose", action="store_true", help="Print errors and skip details")
    args = parser.parse_args()

    zip_path = args.archive_path

    print(f"Counting entries in {zip_path}/peptidelist.txt...")
    expected_count = count_peptidelist_lines(zip_path)
    print(f"Expected entries (lines in index): {expected_count}")

    print(f"\nStarting iterator with all filters set to None (verbose={args.verbose})...")
    actual_count = 0
    for _ in iter_pepbdb(archive_path=zip_path, verbose=args.verbose):
        actual_count += 1
        if actual_count % 1000 == 0:
            print(f"  Yielded {actual_count} entries...")

    print("\n--- Verification Result ---")
    print(f"Expected: {expected_count}")
    print(f"Actual:   {actual_count}")
    print(f"Skipped:  {expected_count - actual_count}")

    if expected_count == actual_count:
        print("\nSUCCESS: Iterator yielded every entry from the index.")
    else:
        print(f"\nWARNING: {expected_count - actual_count} entries were skipped (see verbose output above).")


if __name__ == "__main__":
    main()
