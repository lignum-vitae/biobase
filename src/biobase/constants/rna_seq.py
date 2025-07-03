"""
RNA-Seq Constants and Methods

This module contains constants and methods related to RNA sequencing analysis.

Public module variables:
NORMALIZATION_METHODS -- dictionary containing RNA-seq normalization methods with the following structure:
    'sample' -- methods applied to individual samples (CPM, FPKM, RPKM, TPM)
    'dataset' -- methods applied across entire datasets (DESeq2, edgeR, etc.)
    'batch' -- methods for batch effect correction (ComBat, Harmony, etc.)
Each method contains:
    'method' -- implementation details (R and/or Python packages)
    'description' -- brief description of the method

Key components:
1. Sample-level normalization methods (CPM, FPKM, RPKM, TPM)
2. Dataset-level normalization methods (DESeq2, edgeR, etc.)
3. Batch correction methods (ComBat, Harmony, etc.)

Normalization methods are organized by their scope of application:
- sample: Methods applied to individual samples
- dataset: Methods applied across entire datasets
- batch: Methods specifically for batch effect correction

Each method includes:
- Implementation details (R and/or Python)
- Brief description of the method
- Key references where applicable

Notes:
- Choice of normalization method should be based on experimental design
- Some methods are more appropriate for specific types of analysis
- Batch correction should be applied with caution to avoid removing biological signal
"""

NORMALIZATION_METHODS = {
    "sample": {
        "CPM": {
            "method": "manual",
            "description": "Counts per Million - Scales counts to library size",
        },
        "FPKM": {
            "method": "manual",
            "description": "Fragments per Kilobase of Transcript per Million Fragments Mapped - Accounts for library size and gene length",
        },
        "RPKM": {
            "method": "manual",
            "description": "Reads per Kilobase of Transcript per Million Reads Mapped - Similar to FPKM but for single-end data",
        },
        "TPM": {
            "method": "manual",
            "description": "Transcripts per Million - Similar to FPKM but normalizes at transcript level first",
        },
    },
    "dataset": {
        "size_factors": {
            "method": {"R": "DESeq2", "Python": "PyDESeq2"},
            "description": "Geometric mean normalization accounting for library size differences",
        },
        "TMM": {
            "method": {"R": "edgeR", "Python": "rpy2 + edgeR::calcNormFactors"},
            "description": "Trimmed Mean of M-values - Accounts for RNA composition bias",
        },
        "RLE": {
            "method": {"R": "DESeq2", "Python": "PyDESeq2"},
            "description": "Relative Log Expression - Similar to size factors but uses median instead of mean",
        },
        "quantile": {
            "method": {"R": "limma", "Python": ["scanpy", "sklearn"]},
            "description": "Forces same distribution of counts across samples",
        },
        "upper_quartile": {
            "method": {"R": "edgeR", "Python": "custom"},
            "description": "Normalizes using 75th percentile of counts",
        },
        "median_of_ratios": {
            "method": {"R": "DESeq2", "Python": "PyDESeq2"},
            "description": "DESeq2 default method using median ratios to geometric mean",
        },
        "SCnorm": {
            "method": {"R": "SCnorm", "Python": "rpy2 + SCnorm"},
            "description": "Accounts for systematic bias in count data using quantile regression",
        },
        "scran": {
            "method": {"R": "scran", "Python": "scanpy"},
            "description": "Pool-based size factor estimation for single-cell data",
        },
        "scTransform": {
            "method": {"R": "Seurat", "Python": "rpy2 + Seurat"},
            "description": "Variance stabilizing transformation using regularized negative binomial regression",
        },
        "linnorm": {
            "method": {"R": "linnorm", "Python": "rpy2 + linnorm"},
            "description": "Linear model-based normalization method",
        },
        "VST": {
            "method": {"R": "DESeq2", "Python": "PyDESeq2"},
            "description": "Variance Stabilizing Transformation from DESeq2",
        },
        "rlog": {
            "method": {"R": "DESeq2", "Python": "PyDESeq2"},
            "description": "Regularized log transformation from DESeq2",
        },
        "CSS": {
            "method": {"R": "metagenomeSeq", "Python": "custom"},
            "description": "Cumulative Sum Scaling for sparse count data",
        },
        "UQ": {
            "method": {"R": "edgeR", "Python": "custom"},
            "description": "Upper Quartile normalization",
        },
        "median": {
            "method": "manual",
            "description": "Simple median normalization across samples",
        },
        "PoissonSeq": {
            "method": {"R": "PoissonSeq", "Python": "custom"},
            "description": "Normalization based on Poisson goodness-of-fit statistic",
        },
    },
    "batch": {
        "ComBat": {
            "method": {"R": "sva", "Python": "combat-seq"},
            "description": "Empirical Bayes batch correction",
        },
        "harmony": {
            "method": {"R": "harmony", "Python": "harmonypy"},
            "description": "Fast integration of single-cell data",
        },
        "limma": {
            "method": {"R": "limma", "Python": "custom"},
            "description": "Removes batch effects using linear models",
        },
    },
}
