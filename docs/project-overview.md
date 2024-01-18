# Project Overview

## External dependencies

| Dependency | Use | Installation |
| - | - | - |
| HTSLib     | Reading BAM files.              | [https://www.biostars.org/p/328831/](https://www.biostars.org/p/328831/) <br /> http://www.htslib.org/download/ |
| scikit-learn | Importing the SV scoring model. | [https://www.tensorflow.org/install/lang_c](https://www.tensorflow.org/install/lang_c)|

## SV scoring model

### Datasets
We aim to perform direct comparisons with Sniffles2 results. Here are the
datasets used: [Data
availability](https://www.nature.com/articles/s41587-023-02024-y#data-availability)
 - HG002 PacBio Hifi
    - Saved to data/sniffles2/hifi
 - HG002 ONT (q20_2021.10)
    - Saved to data/sniffles2/ont
 - Medical regions VCF
    - Saved to data/sniffles2/medical_regions
 - 3x 1000 Genomes dipcall benchmark
    - Saved to data/sniffles2/dipcall_benchmark
 - Genome-wide SV v0.6 VCF
    - Saved to data/hg002_giab


### Training the model

## Project structure

```
ContextSV
├── include
│   └── ContextSV.h
├── src
│   └── ContextSV.cpp
└── tests
│   ├── BasicTests.cpp
│   ├── AccuracyTests.cpp
│   └── ScoringModelTests.cpp
└── lib
    ├── tensorflow
    |   ├── LICENSE
    │   ├── include
    │   │   └── ...
    │   └── lib
    │       └── ...
    └── htslib
        ├── LICENSE
        ├── include
        │   └── ...
        └── lib
            └── ...
```
