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
 * [c++ compiler]() It must compile C++11 but C++14 is recommended
 * [gfortran]()
 * [gcc]()
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2

#### Optional dependencies

To have more qp solver:
 * [eigen-qld](https://github.com/jrl-umi3218/eigen-qld.git)
 * [GUROBI](http://www.gurobi.com/) >= 4.0
 * [eigen-gurobi](https://github.com/vsamy/eigen-gurobi)
Also compatible with the LSSOL QP Solver. Unfortunately this not under public license.

To generate the documantation:
 * [Doxygen](http://www.stack.nl/~dimitri/doxygen/)

To have python bindings and unit tests
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.58 (>= 1.21 should work)

#### Building with c++14

```sh
git clone --recursive https://github.com/vsamy/preview_controller
cd preview_controller
./build_and_install
gedit build_and_install_config
./build_and_install
```

Please defines in build_and_install_config where to install the library, the build type, the number of core, etc...
Note that you leave the BOOST_ROOT empty if boost has been installed by default.

##### Specific builds
If you don't want to compile the python bindings you need to
set the variable `PYTHON_BINDINGS` to false (default is true).

If you want to compile the C++ unit tests, you need to set
the variable `BUILD_CXX_TESTS` to true (default is false).

#### Building with c++11

```sh
git clone --recursive https://github.com/vsamy/preview_controller
cd preview_controller
git checkout c++11-back-compatibility
./build_and_install
gedit build_and_install_config
./build_and_install
```

Please defines in build_and_install_config where to install the library, the build type, the number of core, etc...
Note that you leave the BOOST_ROOT empty if boost has been installed by default. 

#### Testing and performance test

You can test the C++ and python version. Those are still basic tests and need to be completed

For c++
```sh
cd _build/tests
./TestSolvers --log_level=all
./TestPreviewControl --log_level=all
```

For python
```sh
cd binding/python/tests
python TestPreviewControl.py
```

## Documentation

You have access to the doxygen files.
Those files are in `<install_path>/share/doc/mpc/doxygen-html/`.
Open the `index.html` file in your web browser.
Plus you can find detailed information [here](https://vsamy.github.io)

## Examples

Please see [here](https://vsamy.github.io) for an example.
You can also check the tests folder.