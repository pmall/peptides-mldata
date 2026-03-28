from peptides_mldata.iterators.pepbdb import iter_pepbdb

def main():
    print("Starting PepBDB exploration and filtering...")

    total_entries = 0
    passed_entries = 0

    iterator = iter_pepbdb(
        archive_path="data/pepbdb-20200318.zip",
        nonstandard_aa=False,
        resolution_min=0,
        resolution_max=2.5,
        mol_type="prot"
    )

    for item in iterator:
        total_entries += 1

        pep_len = item['peptide']['length']
        tgt_len = item['target']['length']

        if not (4 <= pep_len <= 32):
            continue
        if not (15 <= tgt_len <= 128):
            continue

        passed_entries += 1
        print(f"Match: {item['source_key']} | Pep: {pep_len} | Tgt: {tgt_len}")

    print("\n--- Exploration Complete ---")
    print(f"Total entries processed: {total_entries}")
    print(f"Entries passing filters: {passed_entries}")

if __name__ == "__main__":
    main()
