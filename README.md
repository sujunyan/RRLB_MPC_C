# RRLB_MPC_C

This repository contains the generated C code to compare the generated C code from the RRLB MPC controller with some popular MPC solvers. 

## Installation
- CMake: We need cmake to build the code. You can install it using your package manager or from the official website: https://cmake.org/download/
- gcc/g++: We need a C/C++ compiler to build the code. You can install it using your package manager or from the official website: https://gcc.gnu.org/
- HPIPM & BLASFEO: We compare the generated code with HPIPM, which is a high-performance interior point method for quadratic programming. You can find the source code and installation instructions here: https://github.com/giaf/hpipm. 
    - After installing HPIPM, please make sure you run `sudo make install_static` for both HPIPM and BLASFEO to install the static libraries, which are required for linking with our code.
    - Note: If you do not want to run the examples with HPIPM, you can skip this step by commenting out the HPIPM related subdirectory (line 10) in the [CMakeLists.txt](CMakeLists.txt) file.

## How to run the code
1. In the terminal, navigate to the root directory of the project and create a build directory:
```bash
mkdir build
cd build
```

2. Run CMake to configure the project and generate the Makefiles and build files:
```bash
cmake ..
make -j4
```

3. After the build process is complete, you can run the generated executables for each example. We also provide a Python script to run all the examples and compare the results. You can run the script using:
```bash
cd ..   # Navigate back to the root directory
python3 run_examples.py
```

## License
This project is licensed under the BSD 2-Clause License - see the [LICENSE](LICENSE) file for details

Third-party solver license information is documented in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
Local copies/placeholders of third-party license texts are stored in [third_party_licenses/](third_party_licenses/).

## Paper Reference
This codebase is associated with the paper "High-Performance Real-Time Model Predictive Control via Operator Splitting for Embedded Applications" by Yuning Jiang, Junyan Su, Juraj Oravec, and Sorin Olaru. 

