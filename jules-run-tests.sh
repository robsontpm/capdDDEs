#!/bin/bash
set -e

# Setup environment variables needed for build
export BOOST_ROOT=${HOME}/deps
export PATH=${HOME}/deps/bin:$PATH

echo ">> Checking CAPD build..."
if [ ! -f "bin/capd_build/bin/capd-config" ]; then
    echo ">> CAPD not found. Building..."
    if [ -f "./jules-capd-build-nomp.sh" ]; then
        bash ./jules-capd-build-nomp.sh
    else
        echo "Error: jules-capd-build-nomp.sh not found."
    fi
fi

# 2. Rebuild capdDDEs with Coverage enabled
echo ">> Cleaning and configuring build with ENABLE_COVERAGE=ON..."
rm -rf build
mkdir build
cd build

cmake -DENABLE_COVERAGE=ON -DBUILD_PROGRAMS=OFF -DCAPD_DIR=$(pwd)/../bin/capd_build ..

echo ">> Compiling tests..."
make -j4 capd-ddes-DDEForwardTaylorCurvePiece

echo ">> Running tests..."
ctest -V -R capd-ddes-DDEForwardTaylorCurvePiece

echo ">> Generating Coverage..."
# remove mismatch flag for lcov 1.16
lcov --capture --directory . --output-file coverage.info
lcov --extract coverage.info "*/include/capd/ddes/DDEForwardTaylorCurvePiece.h*" --output-file coverage_filtered.info
genhtml coverage_filtered.info --output-directory coverage_report

echo ">> Done."
