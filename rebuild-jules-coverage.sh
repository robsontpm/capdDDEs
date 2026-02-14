#!/bin/bash
# rebuild-jules-coverage.sh
# Script to rebuild the project with coverage enabled, run tests, and generate coverage report.

set -e

WD=$(pwd)
BUILD_DIR="build"
CAPD_DIR="$WD/bin/capd_build"

echo ">> Checking CAPD build..."
if [ ! -f "$CAPD_DIR/bin/capd-config" ]; then
    echo ">> CAPD not found. Building..."
    if [ -f "./jules-capd-build-nomp.sh" ]; then
        bash ./jules-capd-build-nomp.sh
    else
        echo "Error: jules-capd-build-nomp.sh not found. Please ensure CAPD build script is present."
        exit 1
    fi
else
    echo ">> CAPD found in $CAPD_DIR."
fi

# Clean previous build artifacts related to coverage
echo ">> Cleaning build directory..."
if [ -d "$BUILD_DIR" ]; then
    rm -rf "$BUILD_DIR"
fi
mkdir "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with Coverage enabled
echo ">> Configuring CMake with ENABLE_COVERAGE=ON..."
cmake -DENABLE_COVERAGE=ON -DBUILD_PROGRAMS=OFF -DCAPD_DIR="$CAPD_DIR" ..

# Build test targets
echo ">> Compiling tests..."
make -j4 capd-ddes-storage-SharedDoubleton capd-ddes-DDECommon

# Run the tests
echo ">> Running SharedDoubleton tests..."
ctest -V -R "^capd-ddes-storage-SharedDoubleton$"

echo ">> Running DDECommon tests..."
ctest -V -R "^capd-ddes-DDECommon$"

# Generate Coverage Report
echo ">> Generating Coverage Report..."
# Capture coverage data using lcov
# Ignoring mismatch errors which can happen with templates/headers
lcov --capture --directory . --output-file coverage.info --ignore-errors gcov --base-directory ..

# Generate summary
echo ">> Coverage for SharedDoubleton.h:"
lcov --list coverage.info | grep "SharedDoubleton.h"
echo ">> Coverage for DDECommon.h:"
lcov --list coverage.info | grep "DDECommon.h"

echo ">> Done."
