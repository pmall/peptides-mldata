import zipfile
import io
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
    zip_path = "data/pepbdb-20200318.zip"

    print(f"Counting entries in {zip_path}/peptidelist.txt...")
    expected_count = count_peptidelist_lines(zip_path)
    print(f"Expected entries (lines in index): {expected_count}")

    print("\nStarting iterator with all filters set to None (verbose)...")
    actual_count = 0
    for _ in iter_pepbdb(archive_path=zip_path, verbose=True):
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
