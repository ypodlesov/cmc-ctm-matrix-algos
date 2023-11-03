# Basic matrix algorithms

This repository is devoted to study course of computational linear algebra in  **MSU-CMC-CTM** (Moscow State University, Faculty of Computational Mathematics and Cybernetics, Department of Computational Technologies and Modeling).

Here you can also find use-examples of some instruments like **conan package manager**, **google-test**, **google-benchmark** and maybe something else.

## Project file tree
```bash
cmc-ctm-matrix-algos
├── algos
│   ├── benchmarks
│   │   ├── block_vs_parallel.csv
│   │   ├── CMakeLists.txt
│   │   ├── multiply.cpp
│   │   ├── multiply.csv
│   │   ├── parallel_multiply.csv
│   │   └── qr_decomposition.cpp
│   ├── CMakeLists.txt
│   ├── conanfile.txt
│   ├── graphics.ipynb
│   ├── src
│   │   ├── CMakeLists.txt
│   │   ├── common.h
│   │   ├── matrix.cpp
│   │   ├── matrix.h
│   │   ├── qr_decomposition.cpp
│   │   ├── qr_decomposition.h
│   │   ├── vector.cpp
│   │   └── vector.h
│   └── tests
│       ├── CMakeLists.txt
│       ├── main.cpp
│       ├── multiply_base.h
│       ├── multiply_block.cpp
│       ├── multiply_parallel.cpp
│       └── qr_decomposition.cpp
├── .gitignore
└── QR_decomposition_using_boost
    ├── CMakeLists.txt
    ├── conanfile.txt
    └── src
        ├── common.h
        ├── QR.cpp
        ├── QR.h
        └── QR_test.cpp
```

## Structure of main project content components 

1. `algos`
    * `benchmarks` - benchmarking of algorithms
    * `tests` - testing the correctness of algorithm
    * `src` - implementation off all algorithms
2. `QR_decomposition_using_boost` - implementation of **QR decomposition** algorithm using `boost::ublas`.
