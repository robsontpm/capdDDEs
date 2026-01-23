# Building capdDDEs with CMake

This repository now supports building with CMake, alongside the existing Makefiles.

## Prerequisites

*   CMake 3.13 or newer
*   C++17 compliant compiler (GCC, Clang)
*   CAPD library (built and installed, or bundled in `external/`)
*   Boost (optional, for tests)

## Building the Library and Programs

1.  **Build CAPD** (if not already done).
    Follow instructions in `external/README.txt` or use `tldr.sh`.
    Essentially:
    ```bash
    cd external
    ./capd-build.sh
    cd ..
    ```

2.  **Create a build directory** and run CMake:
    ```bash
    mkdir build
    cd build
    cmake ..
    ```

    If CAPD is installed in a non-standard location or you want to use a specific installation, you can hint it:
    ```bash
    cmake .. -DCAPD_DIR=/path/to/capd
    ```
    The build system automatically checks `external/capd` and `bin/capd_build` (the default location for bundled build).

3.  **Compile**:
    ```bash
    cmake --build .
    ```
    or simply `make` (if using Make generator).

    This will build the `capdDDEs` static library and all programs in `programs/examples` and `programs/utils`.

    **Note:** The compiled executables are placed in the `bin` subdirectory of the respective program folder in the source tree (e.g., `programs/examples/mackey-glass-stable-periodic/bin/`).
    To run a program, it is recommended to navigate to that directory first, so that any output files are generated locally:

    ```bash
    cd programs/examples/mackey-glass-stable-periodic/bin
    ./nonrig-find
    ```

## Using capdDDEs in Your Project

A template `CMakeLists.txt` is provided in `docs/CMakeLists.txt`.

To use `capdDDEs` in your own project:

1.  Copy `docs/CMakeLists.txt` to your project directory.
2.  Adjust `CAPDDDES_ROOT` variable in it to point to the `capdDDEs` repository.
3.  Adjust `link_directories` if your build folder is named differently than `build`.
4.  Uncomment and adjust `add_executable` and `target_link_libraries` sections.

Example:
```cmake
add_executable(my_experiment experiment.cpp)
target_link_libraries(my_experiment capdDDEs ${CAPD_LIBRARIES})
```
