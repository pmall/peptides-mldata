from peptides_mldata.iterators.pepbdb import iter_pepbdb

PEPTIDE_MIN_LEN = 4
PEPTIDE_MAX_LEN = 32
SITE_MIN_LEN = 15
SITE_MAX_LEN = 280
RESOLUTION_MIN = 0.001
RESOLUTION_MAX = 2.5
MOL_TYPE = "prot"
NONSTANDARD = False


def main():
    print(f"PepBDB Sweet Spot Count: Pep [{PEPTIDE_MIN_LEN}-{PEPTIDE_MAX_LEN}] | Site [{SITE_MIN_LEN}-{SITE_MAX_LEN}]\n")

    iterator = iter_pepbdb(
        archive_path="data/pepbdb-20200318.zip",
        nonstandard_aa=NONSTANDARD,
        resolution_min=RESOLUTION_MIN,
        resolution_max=RESOLUTION_MAX,
        mol_type=MOL_TYPE,
    )

    total = 0
    passed = 0
    for item in iterator:
        total += 1
        peptide_len = item['peptide']['length']
        site_len = item['target']['length']
        if PEPTIDE_MIN_LEN <= peptide_len <= PEPTIDE_MAX_LEN and SITE_MIN_LEN <= site_len <= SITE_MAX_LEN:
            passed += 1
        if total % 1000 == 0:
            print(f"  {total} entries scanned, {passed} matches...")

    print(f"\nTotal scanned: {total}")
    print(f"Sweet spot entries: {passed} ({passed/total*100:.1f}%)")


if __name__ == "__main__":
    main()
