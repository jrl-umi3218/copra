# Copra

Copra (**Co**ntrol & **pr**eview **a**lgorithms) is a C++ library implementing
linear model predictive control. It relies on quadratic programming (QP)
solvers. Python bindings are available.

## Installation

Compilation has been tested on Linux (gcc/clang/msvc).

### Dependencies

* Any compiler with C++14 support
* [CMake](https://cmake.org) >= 2.8
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/): to generate documentation
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
* [Git](https://git-scm.com/)
* [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
* [eigen-quadprog](https://github.com/vsamy/eigen-quadprog)

#### Optional

* [Catch2](https://github.com/catchorg/Catch2): for tests
* [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.58: for Python bindings and unit tests
* [GUROBI](http://www.gurobi.com/) >= 4.0: optional QP solver
* [eigen-gurobi](https://github.com/vsamy/eigen-gurobi): bindings for GUROBI
* [eigen-qld](https://github.com/jrl-umi3218/eigen-qld.git): optional QP solver
* [eigen-osqp](https://github.com/jrl-umi3218/eigen-osqp.git): optional QP solver
* [pygen-converter](https://github.com/vsamy/pygen-converter): for python bindings

### Building from source on Linux

The library is written in c++14. Compiler that can not support it won't be able to compile it.

```sh
git clone --recursive git@github.com:vsamy/copra.git
cd Copra
mkdir build && cd build
cmake ..
ccmake .  # configure e.g. PYTHON_BINDINGS or BUILD_CXX_TESTS
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
``BUILD_CXX_TESTS`` option:

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

Doxygen files will be compiled into
``<install_path>/share/doc/mpc/doxygen-html/``. Open the `index.html` file in
your web browser. You can also check out unit tests, as well as the following
two examples:

* [C++ example of Copra](https://vsamy.github.io/en/blog/copra-example-cpp)
* [Python example of Copra](https://vsamy.github.io/en/blog/copra-example-python)
