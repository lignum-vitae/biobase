# Biobase

[![Static Badge](https://img.shields.io/badge/Project_Name-Biobase-blue)](https://github.com/lignum-vitae/biobase)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

A Python package providing standardized biological constants and scoring matrices for bioinformatics pipelines. Biobase aims to eliminate the need to repeatedly recreate common biological data structures and scoring systems in your code.

## Table of Contents
- [Quick Start](#quick-start)
  - [Access amino acid properties](#access-amino-acid-properties)
  - [Use scoring matrices](#use-scoring-matrices)
  - [Analyze DNA sequences](#analyze-dna-sequences)
  - [Find protein motifs](#find-protein-motifs)
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

#### Access amino acid properties:
```python
from biobase.constants import ONE_LETTER_CODES, MONO_MASS
print(ONE_LETTER_CODES)  # 'ACDEFGHIKLMNPQRSTVWY'
print(MONO_MASS['A'])    # 71.037113805
```
#### Use scoring matrices:
```python
from biobase.matrix import Blosum
blosum62 = Blosum(62)
print(blosum62['A']['A'])  # 4
print(blosum62['W']['C'])  # -2
```
#### Analyze DNA sequences:
```python
from biobase.analysis import Dna
sequence = "ATCGTAGC"
print(Dna.transcribe(sequence))         # 'AUCGUAGC'
print(Dna.complement_dna(sequence))     # 'GCTACGAT'
print(Dna.calculate_gc_content(sequence))  # 50.0
```
#### Find protein motifs:
```python
from biobase.analysis import find_motifs
sequence = "ACDEFGHIKLMNPQRSTVWY"
print(find_motifs(sequence, "DEF"))  # [3]
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

Biobase aims to provide Python-friendly versions of common biological constants and tools for bioinformatics pipelines. Key objectives:

1. Standardize biological data structures
2. Provide efficient implementations of common scoring systems
3. Ensure type safety and validation
4. Maintain comprehensive documentation
5. Support modern Python practices

## Contributing

We welcome contributions! Please read our:
- [Code of Conduct](https://github.com/lignum-vitae/biobase/blob/main/docs/CODE_OF_CONDUCT.md)
- [Contribution Guidelines](https://github.com/lignum-vitae/biobase/blob/main/docs/CONTRIBUTING.md)

## Project Status

### Current Version: 0.4.1-alpha

#### Core Features
- ✅ BLOSUM and PAM matrix implementations
- ✅ Basic amino acid constants and conversions
- ✅ DNA/RNA sequence analysis tools
- ✅ Protein motif searching
- ✅ Core biological constants
- ✅ Additional scoring matrices
- ✅ Extended amino acid properties
  
#### Analysis Tools
- ✅ GC content calculation
- ✅ DNA/RNA transcription
- ✅ DNA complementation
- ✅ Motif finding
- 🚧 File format parsers (FASTA, GenBank, etc.)
- 📋 Statistical analysis tools


#### Documentation
- ✅ Basic README
- ✅ Code of Conduct
- ✅ Contributing Guidelines
- ✅ Usage Examples

#### Development
- 🚧 PyPI package deployment
- 🚧 CI/CD Pipeline
- 🚧 Code Coverage
- 📋 Automated Releases

### Legend
- ✅ Complete
- 🚧 In Progress
- 📋 Planned

### Stability
This project is in the alpha stage. APIs may change without warning until version 1.0.0.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/lignum-vitae/biobase/blob/main/LICENSE) file for details.
