import argparse
import csv
from peptides_mldata.iterators.pepbdb import iter_pepbdb

# Iterator quality defaults
RESOLUTION_MIN = 0.001
RESOLUTION_MAX = 2.5
MOL_TYPE = "prot"
NONSTANDARD = False


def pepbdb_metadata(min_pep_len: int, max_pep_len: int, min_target_len: int, max_target_len: int, output_csv: str):
    print(f"Exporting PepBDB metadata to {output_csv}...")
    print(f"Filters: Pep [{min_pep_len}-{max_pep_len}] | Site [{min_target_len}-{max_target_len}]\n")

    iterator = iter_pepbdb(
        archive_path="data/pepbdb-20200318.zip",
        nonstandard_aa=NONSTANDARD,
        resolution_min=RESOLUTION_MIN,
        resolution_max=RESOLUTION_MAX,
        mol_type=MOL_TYPE,
        min_pep_len=min_pep_len,
        max_pep_len=max_pep_len,
        min_target_len=min_target_len,
        max_target_len=max_target_len,
    )

    count = 0
    with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
        fieldnames = ["pdb_id", "peptide_chain", "peptide_len", "target_chain", "target_len", "source_key"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for item in iterator:
            writer.writerow({
                "pdb_id": item["pdb_id"],
                "peptide_chain": item["peptide"]["chain"],
                "peptide_len": item["peptide"]["length"],
                "target_chain": item["target"]["chain"],
                "target_len": item["target"]["length"],
                "source_key": item["source_key"]
            })
            count += 1
            if count % 1000 == 0:
                print(f"  {count} entries exported...")

    print(f"\nExport complete. Total entries: {count}")


def main():
    parser = argparse.ArgumentParser(description="Export filtered PepBDB metadata to a CSV file.")
    parser.add_argument("output", help="Path to the output CSV file")
    parser.add_argument("--min_pep_len", type=int, default=4, help="Minimum peptide length (default: 4)")
    parser.add_argument("--max_pep_len", type=int, default=32, help="Maximum peptide length (default: 32)")
    parser.add_argument("--min_target_len", type=int, default=15, help="Minimum target/site length (default: 15)")
    parser.add_argument("--max_target_len", type=int, default=280, help="Maximum target/site length (default: 280)")
    args = parser.parse_args()

    pepbdb_metadata(
        args.min_pep_len,
        args.max_pep_len,
        args.min_target_len,
        args.max_target_len,
        args.output
    )


if __name__ == "__main__":
    main()
