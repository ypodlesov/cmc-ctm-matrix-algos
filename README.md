# Basic matrix algorithms

This repository is devoted to study course of computational linear algebra in  **MSU-CMC-CTM** (Moscow State University, Faculty of Computational Mathematics and Cybernetics, Department of Computational Technologies and Modeling).

Here you can also find use-examples of some instruments like **conan package manager**, **google-test**, **google-benchmark** and maybe something else.

## Project file tree

`git ls-tree -r --name-only HEAD | tree --fromfile`

```bash
cmc-ctm-matrix-algos
├── CMakeLists.txt
├── .gitignore
├── qr_decomposition
│   ├── CMakeLists.txt
│   ├── common.h
│   ├── matrix.cpp
│   ├── matrix.h
│   ├── qr_decomposition.cpp
│   ├── qr_decomposition.h
│   ├── run.cpp
│   ├── test.cpp
│   ├── vector.cpp
│   └── vector.h
├── README.md
└── tools
    ├── catch_main.cpp
    ├── cmake
    │   ├── BuildFlags.cmake
    │   ├── FindCatch.cmake
    │   ├── FindFFTW.cmake
    │   ├── Protobuf.cmake
    │   └── TestSolution.cmake
    └── commons
        ├── run_channel.h
        ├── runner.h
        ├── test_channel.h
        └── util.h
```