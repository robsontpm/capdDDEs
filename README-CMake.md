# Building capdDDEs with CMake

This repository now supports building with CMake, alongside the existing Makefiles.

This support is now experimental and might change in the near future.

In a long run cmake should replace old method (Makefiles).

## Prerequisites

* CMake 3.13 or newer
* C++17 compliant compiler (GCC, Clang)
* CAPD library (built and installed, or bundled in `external/`)
* Boost (optional, for tests)

## Building the Library and Programs

1. **Create a build directory** and run CMake:

    ```bash
    mkdir build
    cd build
    cmake ..
    ```

    In the configuration step, you can optionally pass the following options:
    * `-DcapdDDEs_BUILD_TESTS=<ON|OFF>`  
      enable/disable building tests. (default: `ON`)
    * `-DcapdDDEs_BUILD_EXAMPLES=<ON|OFF>`  
      enable/disable building examples. (default: `ON`)
    * `-DcapdDDEs_BUILD_UTILS=<ON|OFF>`  
      enable/disable building utility programs. (default: `ON`)
    * `-DcapdDDEs_USE_SYSTEM_CAPD=<ON|OFF>`  
      enable/disable searching for the system installation of CAPD library.
      (default: `ON`)
    * `-DcapdDDEs_INSTALL=<ON|OFF>`
      enable/disable installation of the library when calling `cmake --install`
      potentially useful when `capdDDEs` is meant to be used as a private
      dependency when added through `add_subdirectory` cmake command.
      If enabled, the files will be installed to the path specified by
      `CMAKE_INSTALL_PREFIX`. (default: `ON`)

    When using system `CAPD`, you can additionally specify the path to the
    `find-capd.cmake` file through the `-Dcapd_DIR=<path/to/capd>` option.

    For example, to build `capdDDEs` library with `CAPD` installed in `~/.local`
    without building examples and utility programs, but building tests, you
    could use the following command:

    ```bash
    cmake ..\
      -DcapdDDEs_BUILD_TESTS=ON\
      -DcapdDDEs_BUILD_EXAMPLES=OFF\
      -DcapdDDEs_BUILD_UTILS=OFF\
      -DcapdDDEs_USE_SYSTEM_CAPD=OFF\
      -Dcapd_DIR=~/.local/lib/cmake/capd
    ```

2. **Compile**:

    ```bash
    cmake --build .
    ```

    or simply `make` (if using Make generator).

    This will build the `capdDDEs` static library and all enabled additional
    components

    **Note:** The compiled executables are placed in the `bin` subdirectory of
    the respective program folder in the source tree (e.g.,
    `programs/examples/mackey-glass-stable-periodic/bin/`).
    To run a program, it is recommended to navigate to that directory first, so
    that any output files are generated locally:

    ```bash
    cd programs/examples/mackey-glass-stable-periodic/bin
    ./nonrig-find
    ```

3. **Install** (*optional*)

    ```bash
    cmake --install .
    ```

    This command will install the library to the path specified by the
    `CMAKE_INSTALL_PREFIX` variable set during the configuration phase.

## Using capdDDEs in Your Project

A template `CMakeLists.txt` is provided in `docs/CMakeLists.txt`.

To use `capdDDEs` in your own project you can either:

* If `capdDDEs` is already installed, call `find_package(capdDDEs)` (it might
  be necessary to set `capdDDEs_DIR` variable to a directory containing
  `capdddes-config.cmake` file).
* Clone the `capdDDEs` repository into your project and call
  `add_subdirectory(capdDDEs)`.
* Use the cmake `FetchContent` module.

All of these options define a `capdDDEs::capdDDEs` target. You can link with
`capdDDEs` by calling `target_link_libraries(<target> ... capdDDEs::capdDDEs)`.

Example:

```cmake
add_executable(my_experiment experiment.cpp)
target_link_libraries(my_experiment PRIVATE capdDDEs::capdDDEs)
```
