# Copra

Copra (**Co**ntrol & **pr**eview **a**lgorithms) is a C++ library implementing
linear model predictive control. It relies on quadratic programming (QP)
solvers. Python bindings are available.

## Installing

### Dependencies

* [CMake]() >= 2.8
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
* [Git]()
* [c++ compiler]() It must compile C++11 but C++14 is recommended
* [doxygen]()
* [gcc]()
* [gfortran]()
* [pkg-config]()

#### Optional

* [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.58, for Python bindings and unit tests
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) to generate documentation
* [GUROBI](http://www.gurobi.com/) >= 4.0, optional QP solver
* [eigen-gurobi](https://github.com/vsamy/eigen-gurobi), bindings for GUROBI
* [eigen-qld](https://github.com/jrl-umi3218/eigen-qld.git), optional QP solver

### Building

The library assumes you are compiling C++14 by default:

```sh
git clone --recursive git@github.com:stephane-caron/Copra.git
cd Copra
mkdir build && cd build
cmake ..
ccmake .  # configure e.g. PYTHON_BINDINGS or BUILD_CXX_TESTS
make -j4
sudo make install
```

Take a look at the ``c++11-back-compatibility`` branch if you'd rather compile
for C++11.

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

- [C++ example of Copra](https://vsamy.github.io/en/blog/copra-example-cpp)
- [Python example of Copra](https://vsamy.github.io/en/blog/copra-example-python)
