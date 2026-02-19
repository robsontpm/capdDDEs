#!/bin/bash
# rebuild-jules-tests.sh
# Script to reproduce the testing workflow for capdDDEs, including environment setup and test execution.
# This workflow reproduces the std::bad_alloc error in SharedDoubleton tests when they are enabled.

set -e

# 1. Build CAPD (if not already present in bin/capd_build)
# Using jules-capd-build-nomp.sh as per previous steps (MPFR off) or default if preferred.
# For exact reproduction of my environment, I used jules-capd-build-nomp.sh.
echo ">> Checking CAPD build..."
if [ ! -f "bin/capd_build/bin/capd-config" ]; then
    echo ">> CAPD not found. Building..."
    # Ensure jules-capd-build-nomp.sh exists (created in previous steps)
    if [ -f "./jules-capd-build-nomp.sh" ]; then
        bash ./jules-capd-build-nomp.sh
    else
        echo "Error: jules-capd-build-nomp.sh not found. Please ensure CAPD build script is present."
        exit 1
    fi
else
    echo ">> CAPD found in bin/capd_build."
fi

# 2. Rebuild capdDDEs with Coverage enabled
echo ">> Cleaning and configuring build with ENABLE_COVERAGE=ON..."
rm -rf build
mkdir build
cd build

# Compiler flags observed during failure:
# -std=gnu++17 -std=c++17 -g -O0 --coverage -frounding-math ...
# These are set by cmake -DENABLE_COVERAGE=ON and the project configuration.

export BOOST_ROOT="${HOME}/deps/boost"
cmake -DENABLE_COVERAGE=ON -DBUILD_PROGRAMS=OFF -DCAPD_DIR=$(pwd)/../bin/capd_build -DBOOST_ROOT="$BOOST_ROOT" ..

echo ">> Compiling..."
make -j4

# 3. Run all tests
ctest -V -R 

echo ">> Done."
