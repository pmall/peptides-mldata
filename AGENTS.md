# Project peptides-mldata

## Tech Stack

- **Language**: Python
- **Package Manager**: uv (MANDATORY: Do NOT use `pip` directly)

## Goal

Extract clean peptide sequence and binding site pairs from various sources (mostly processed PDB data) to serve as a dataset for training machine learning models in downstream projects. The data will be consolidated into an LMDB file for efficient access.

## Data Extraction Rules

- **Extraction targets**:
  - Amino acid sequences
  - 3D coordinates, restricted strictly to the **3 backbone atoms** (N, CA, C) per residue instead of the full 37.
- **Workflow**:
  - Store source data as compressed files (tar/zip) in the ignored `data/` directory.
  - Parse data using **on-the-fly iterators** to avoid extracting source datasets and cluttering the filesystem.
  - Iterators are responsible for yielding a unique, source-specific key per entry (e.g. `1a1r_C` for PepBDB).
  - **Iterator Parameter Philosophy**:
    - Expose metadata (e.g., resolution, mol_type) and filtering constraints (e.g., `min_pep_len`, `max_pep_len`, `min_target_len`, `max_target_len`) as optional parameters (defaulting to `None`) on the iterator function itself.
    - **Perform filtering as early as possible**: If a dataset provides a metadata index (e.g., a text file inside the zip), use it to skip entries before parsing expensive raw data (like PDB coordinates) to optimize performance.
- **Script Architecture Rules**:
  - **Functional Separation**: Separate CLI argument parsing from execution logic. Extract the core processing logic into a standalone function (e.g., `process_datasource(min_pep_len, ...)`) that takes the parameters directly.
  - **Argparse Integration**: Use `argparse` in `main()` to handle user inputs, providing sensible defaults (e.g., peptide length 4-32) for standard "sweet spot" processing.
- **Nullability and Strictness**:
  - Prefer strong typing and exact data representation.
  - **Error Handling in Iterators**: Since Python generators cannot resume after raising an exception, iterators must catch exceptions internally and skip malformed entries to keep iterating. All iterators must expose a `verbose: bool = False` parameter. When `verbose` is `True`, exceptions are printed to stdout. When `False`, entries are skipped silently. This keeps the iterator alive while giving the caller visibility into data quality issues.
  - Do NOT preemptively set every field to nullable. Start strict, and adjust the schema structure as null values are discovered during parsing.

## LMDB Strategy

- **One LMDB per source**: each `*_build.py` script writes to a dedicated LMDB (e.g. `data/pepbdb.lmdb`). Keys are the bare `source_key` returned by the iterator — **no source prefix is added**.
- **Future merge step**: a separate script will be responsible for deduplicating and merging per-source LMDBs into a unified dataset, handling key conflicts across sources at that point.

## Output Schema

The data will be processed and stored in an LMDB file using the following core schema format:

```json
{
    "source_db": "pepbdb",
    "pdb_id": "1A2B",
    "target": {
        "chain": "A",
        "length": N,                     // Integer length of binding site
        "taxon_ids": ["9606"],           // List of NCBI Taxonomy IDs
        "protein_ids": ["P12345"],       // List of UniProt accessions
        "sequence": "MVLSPADK...",       // binding site of length N
        "3d_coordinates": [N, 3]         // Numpy float32 array
    },
    "peptide": {
        "chain": "C",
        "length": P,                     // Integer length of peptide
        "taxon_ids": ["111938"],         // List of NCBI Taxonomy IDs
        "protein_ids": ["P98765"],       // List of UniProt accessions
        "sequence": "YGGFL",             // Peptide sequence of length P
        "3d_coordinates": [P, 3]         // Numpy float32 array
    }
}
```
