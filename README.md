# preview_controller

preview_controller is a library implementing a model preview controller.
It uses QP softwares to solve its problem. A python binding is available.

## Installing

### Manual

#### Dependencies

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 2.8
 * [pkg-config]()
 * [doxygen]()
 * [c++ compiler]() Version to compile C++14 (g++ => 5.0, clang++ >= 3.4)
 * [gfortran]()
 * [gcc]()
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.58 (>= 1.49 may work)
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2

#### Optional dependencies

To have more qp solver:
 * [GUROBI](http://www.gurobi.com/) >= 4.0
 * [eigen-gurobi](https://github.com/jrl-umi3218/eigen-gurobi)

#### Building

```sh
git clone --recursive https://github.com/vsamy/preview_controller
cd preview_controller
./build_and_install
gedit build_and_install_config
./build_and_install
```

Please defines in build_and_install_config where to install the library, the build type, the number of core, etc...
Note that you leave the BOOST_ROOT empty if boost has been installed by default. 

Where the main options are:

#### Testing and performace test

You can test the C++ and python version. Those are still basic tests and need to be completed

For c++
```sh
cd _build/tests
./TestSolvers
./TestPreviewControl --log_level=all
```

For python
```sh
cd binding/python/tests
python TestPreviewControl.py
```

#### Examples

There is not basic examples yet. Please see test files for an overview.
Please see the doxygen files for the documentation.
