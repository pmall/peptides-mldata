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
  - Store source data as compressed files (zip) in the ignored `data/` directory.
  - Parse data using **on-the-fly iterators** to avoid extracting source datasets and cluttering the filesystem.
  - Iterators are responsible for yielding a unique, source-specific key per entry. The final LMDB storage keys will be composite (e.g., `source:source_key`) to prevent collisions across datasets.
  - **Iterator Parameter Philosophy**: Do not muddy the standardized output schema with dataset-specific metadata. If a dataset source has unique metadata (like resolution or mol_type), expose them exclusively as optional parameters (defaulting to `None`) on the iterator function itself to allow upstream filtering!
- **Nullability and Strictness**:
  - Prefer strong typing and exact data representation over defensive programming.
  - Do NOT preemptively set every field to nullable. Start strict, and adjust the schema structure as null values are discovered during parsing.

## Output Schema

The data will be processed and stored in an LMDB file using the following core schema format:

```json
{
    "source_db": "pepbdb",
    "pdb_id": "1A2B",
    "target": {
        "chain": "A",
        "length": N,                     // Integer length of binding site
        "taxon_id": null,                // e.g., 9606 for Human, or null if missing in source
        "protein_id": null,              // e.g., UniProt ID P12345, or null if missing in source
        "sequence": "MVLSPADK...",       // binding site of length N
        "3d_coordinates": [N, 3]         // Numpy float32 array
    },
    "peptide": {
        "chain": "C",
        "length": P,                     // Integer length of peptide
        "taxon_id": null,                // e.g., 111938 for Viral, or null if synthetic/missing
        "protein_id": null,              // e.g., Parent UniProt ID P98765, or null if missing
        "sequence": "YGGFL",             // Peptide sequence of length P
        "3d_coordinates": [P, 3]         // Numpy float32 array
    }
}
```
