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
  - Iterators are responsible for yielding a unique, source-specific key per entry, avoiding a strict global key convention.
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
        "taxon_id": 9606,                // e.g., Human
        "protein_id": "P12345",          // UniProt ID
        "sequence": "MVLSPADK...",       // binding site of length N
        "3d_coordinates": [N, 3]         // Numpy float32 array
    },
    "peptide": {
        "chain": "C",
        "taxon_id": 111938,              // e.g., Viral (or null if synthetic)
        "protein_id": "P98765",          // Parent UniProt ID (often null)
        "sequence": "YGGFL",             // Peptide sequence of length P
        "3d_coordinates": [P, 3]         // Numpy float32 array
    }
}
```
