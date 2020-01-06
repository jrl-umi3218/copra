# Copra

[![License](https://img.shields.io/badge/License-BSD%202--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![CI](https://github.com/jrl-umi3218/copra/workflows/CI%20of%20copra/badge.svg?branch=master)](https://github.com/jrl-umi3218/copra/actions?query=workflow%3A%22CI+of+copra%22)
[![Documentation](https://img.shields.io/badge/doxygen-online-brightgreen?logo=read-the-docs&style=flat)](http://jrl-umi3218.github.io/copra/doxygen/HEAD/index.html)

Copra (**Co**ntrol & **pr**eview **a**lgorithms) is a C++ library implementing
linear model predictive control. It relies on quadratic programming (QP)
solvers. Python bindings are available.

This work was originally made by [Vincent Samy](https://github.com/vsamy).

Copra is licensed under the [BSD-2-Clause](https://opensource.org/licenses/BSD-2-Clause). However, its default QP solver (eigen-quadprog) is licensed under the LGPL-2 and this cannot be changed. Please be aware of the related restrictions if you plan to work with copra. At a later date, we might switch copra default QP solver to one with a less restrictive license.

## Installation

Copra should compiles and is tested on Linux, macOS and Windows.

### Ubuntu LTS (16.04, 18.04, 20.04)

```bash
# Make sure you have required tools
sudo apt install apt-transport-https lsb-release ca-certificates gnupg
# Add our key
sudo apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key 892EA6EE273707C6495A6FB6220D644C64666806
# Add our repository (stable versions)
sudo sh -c 'echo "deb https://dl.bintray.com/gergondet/multi-contact-release $(lsb_release -sc) main" | sudo tee /etc/apt/sources.list.d/multi-contact.list'
# Use this to setup the HEAD version
# sudo sh -c 'echo "deb https://dl.bintray.com/gergondet/multi-contact-head $(lsb_release -sc) main" | sudo tee /etc/apt/sources.list.d/multi-contact.list'
# Update packages list
sudo apt update
# Install packages
sudo apt install libcopra-dev
```

### Dependencies

* C++ compiler with C++14 support (see below for C++11 support)
* [CMake](https://cmake.org) >= 2.8
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/): to generate documentation (optional)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
* [eigen-quadprog](https://github.com/vsamy/eigen-quadprog)

#### Optional

* [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) for Python bindings (>= 1.58) and unit tests
* [GUROBI](http://www.gurobi.com/) >= 4.0: optional QP solver
* [eigen-gurobi](https://github.com/vsamy/eigen-gurobi): bindings for GUROBI
* [eigen-qld](https://github.com/jrl-umi3218/eigen-qld.git): optional QP solver
* [eigen-osqp](https://github.com/jrl-umi3218/eigen-osqp.git): optional QP solver
* [pygen-converter](https://github.com/vsamy/pygen-converter): for python bindings

### Building from source on Linux

The library is written in c++14. Compiler that can not support it won't be able to compile it.

```sh
git clone --recursive git@github.com:jrl-umi3218/copra.git
cd Copra
mkdir build && cd build
cmake ..
ccmake .  # configure e.g. PYTHON_BINDING or BUILD_TESTING
make -j4
sudo make install
```

### C+11 compatible code

A c++11 compatible is available on a specific branch.
Although it exists, it may not be updated as often as the master branch.
The branch is called `c++11-compatible`. To use it, before compilation and after cloning the repo, do:

```sh
git checkout c++11-compatible
```

Then follow the steps in the section just above.

### Testing

C++ tests will be compiled in your build folder if you enabled the
``BUILD_TESTING`` option:

```sh
cd build/tests
./TestSolvers --log_level=all
./TestLMPC --log_level=all
```

Once Python bindings are installed, you can check them with:

```sh
cd binding/python/tests
python TestMPController.py
```

## Documentation

Doxygen documentation is available [online](http://jrl-umi3218.github.io/copra/doxygen/HEAD/index.html)

You can also check out unit tests, as well as the following two examples:

* [C++ example of Copra](https://vsamy.github.io/en/blog/copra-example-cpp)
* [Python example of Copra](https://vsamy.github.io/en/blog/copra-example-python)
