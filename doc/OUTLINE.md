# ContextSV Outline

## External dependencies

| Dependency | Use | Installation |
| - | - | - |
| HTSLib     | Reading BAM files.              | [https://www.biostars.org/p/328831/](https://www.biostars.org/p/328831/) |
| Tensorflow | Importing the SV scoring model. | |

## SV scoring model
The SV scoring model is trained and saved as a Tensorflow weighted graph in PB file format. It can be trained via the [Python API](https://www.tensorflow.org/api_docs/python/tf). The [C++ API](https://www.tensorflow.org/api_docs/cc) is used to import and utilize this model for scoring SVs.

### Training the model

## Project structure

```
ContextSV
├── include
│   └── ContextSV
│       └── Caller.h
├── src
│   ├── ContextSV.cpp
│   ├── Integrator.h
│   └── Integrator.cpp
└── tests
│   ├── CallerTests.cpp
│   └── IntegratorTests.cpp
└── lib
    ├── tensorflow
    │   ├── include
    │   │   └── ...
    │   ├── src
    │   │   └── ...
    │   └── test
    │       └── ...
    └── htslib
        ├── include
        │   └── sam.h
```

