# ContextSV Outline

## External dependencies

| Dependency | Use | Installation |
| - | - | - |
| HTSLib     | Reading BAM files.              | [https://www.biostars.org/p/328831/](https://www.biostars.org/p/328831/) <br /> http://www.htslib.org/download/ |
| Tensorflow | Importing the SV scoring model. | [https://www.tensorflow.org/install/lang_c](https://www.tensorflow.org/install/lang_c)|

## SV scoring model
The SV scoring model is trained and saved as a Tensorflow weighted graph in PB file format. It can be trained via the [Python API](https://www.tensorflow.org/api_docs/python/tf). The [C++ API](https://www.tensorflow.org/api_docs/cc) is used to import and utilize this model for scoring SVs.

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

## Tasks by order
- [x] Draft project structure
- [ ] CLI with BAM file input
