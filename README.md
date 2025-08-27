# Biobase

[![Static Badge](https://img.shields.io/badge/Project_Name-Biobase-blue)](https://github.com/lignum-vitae/biobase)
[![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Flignum-vitae%2Fbiobase%2Fmain%2Fpyproject.toml)](https://github.com/lignum-vitae/biobase/blob/main/pyproject.toml)
[![PyPI version](https://img.shields.io/pypi/v/biobase.svg)](https://pypi.python.org/pypi/biobase)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![GitHub branch status](https://img.shields.io/github/checks-status/lignum-vitae/biobase/main)](https://github.com/lignum-vitae/biobase)

A Python package providing standardized biological constants and scoring matrices
for bioinformatics pipelines.
Biobase aims to eliminate the need to repeatedly recreate common biological data
structures and scoring systems in your code.

## Table of Contents

- [Quick Start](#quick-start)
  - [Access amino acid properties](#access-amino-acid-properties)
  - [Use scoring matrices](#use-scoring-matrices)
  - [Analyze DNA sequences](#analyze-dna-sequences)
  - [Find protein motifs](#find-protein-motifs)
  - [Parse FASTA](#parse-fasta)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Regular Installation](#regular-installation)
  - [Development Installation](#development-installation)
- [Running Files](#running-files)
- [Data Files](#data-files)
- [Project Goals](#project-goals)
- [Contributing](#contributing)
- [Project Status](#project-status)
- [License](#license)

## Quick Start

#### Access amino acid properties

```python
from biobase.constants import ONE_LETTER_CODES, MONO_MASS
print(ONE_LETTER_CODES)  # 'ACDEFGHIKLMNPQRSTVWY'
print(MONO_MASS['A'])    # 71.037113805
`FastaRecord``
```

#### Use scoring matrices

```python
from biobase.matrix import Blosum
blosum62 = Blosum(62)
print(blosum62['A']['A'])  # 4
print(blosum62['W']['C'])  # -2
```

#### Analyze DNA sequences

```python
from biobase.analysis import Dna
sequence = "ATCGTAGC"
print(Dna.transcribe(sequence))               # 'AUCGUAGC'
print(Dna.complement(sequence))               # 'TAGCATCG'
print(Dna.complement(sequence, reverse=True)) # 'GCTACGAT'
print(Dna.calculate_gc_content(sequence))     # 50.0
print(Dna.entropy(sequence))                  # 2.0

seq = "ccatgccctaaatggggtag"
for start, end, orf in Dna.find_orfs(seq, include_seq=True)
    print(start, end, orf)
# 2, 11, "ATGCCCTAA"
# 11, 20, "ATGGGGTAG"
```

#### Find protein motifs

```python
from biobase.analysis import find_motifs
sequence = "ACDEFGHIKLMNPQRSTVWY"
print(find_motifs(sequence, "DEF"))  # [3]
```

#### Parse FASTA

```python
from biobase.parser import fasta_parser
records = fasta_parser(fasta)
for r in records:
    print(r.id) # CAA39742.1
    print(r.seq) # MTNIRKSHPLMKII...
```

## Requirements

- Python 3.10+
- pip (for installation)

## Installation

### Regular Installation

`pip install biobase`

### Development Installation

Clone the repository and install in editable mode:

```nginx
git clone https://github.com/lignum-vitae/biobase.git
cd biobase
pip install -e .
```

## Running Files

To ensure relative imports work correctly, always run files using the module path from the project root:

### Run a specific file

python -m src.biobase.matrix

## Data Files

- `src/biobase/matrices/`: Scoring matrix data stored in JSON file format

## Project Goals

Biobase aims to provide Python-friendly versions of common biological constants
and tools for bioinformatics pipelines. Key objectives:

1. Standardize biological data structures
2. Provide efficient implementations of common scoring systems
3. Ensure type safety and validation
4. Maintain comprehensive documentation
5. Support modern Python practices

## Contributing

We welcome contributions! Please read our:
- [Code of Conduct](https://github.com/lignum-vitae/biobase/blob/main/docs/CODE_OF_CONDUCT.md)
- [Contribution Guidelines](https://github.com/lignum-vitae/biobase/blob/main/docs/CONTRIBUTING.md)

### Stability

This project is in the beta stage. APIs may change without warning until version 1.0.0.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/lignum-vitae/biobase/blob/main/LICENSE) file for details.
