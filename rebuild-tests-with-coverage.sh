#!/bin/bash
# rebuild-jules-coverage.sh
# Script to rebuild the project with coverage enabled, run SharedDoubleton tests, and generate coverage report.
# Similar to rebuild-jules-tests.sh but focused on coverage analysis.

set -e

WD=$(pwd)
BUILD_DIR="build"
CAPD_DIR="$WD/bin/capd_build"

# Define lcov command
if [ -f "$WD/_deps/lcov-1.16/bin/lcov" ]; then
    LCOV_CMD="$WD/_deps/lcov-1.16/bin/lcov"
    echo ">> Using local lcov: $LCOV_CMD"
else
    LCOV_CMD="lcov"
    echo ">> Using system lcov: $LCOV_CMD"
fi

# Define Boost options if local boost exists
BOOST_OPTS=""
if [ -d "$WD/_deps/boost_1_82_0" ]; then
    echo ">> Using local Boost at $WD/_deps/boost_1_82_0"
    BOOST_OPTS="-DBOOST_ROOT=$WD/_deps/boost_1_82_0 -DBoost_ROOT=$WD/_deps/boost_1_82_0 -DBoost_INCLUDE_DIR=$WD/_deps/boost_1_82_0 -DBoost_NO_SYSTEM_PATHS=ON"
fi

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
cmake -DENABLE_COVERAGE=ON -DBUILD_PROGRAMS=OFF -DCAPD_DIR="$CAPD_DIR" $BOOST_OPTS ..

# make all
make -j4

# Run the tests
echo ">> Running ALL tests..."
ctest -V

# Generate Coverage Report
echo ">> Generating Coverage Report..."
# Capture coverage data using lcov
# Ignoring mismatch errors which can happen with templates/headers
$LCOV_CMD --capture --directory . --output-file coverage.info --ignore-errors gcov --base-directory ..

# Generate summary
echo ">> Coverage for specific components:"
$LCOV_CMD --list coverage.info | grep "DDECommon.h"
$LCOV_CMD --list coverage.info | grep "DDECommon.cpp"
$LCOV_CMD --list coverage.info | grep "DoubletonInterface.h"
$LCOV_CMD --list coverage.info | grep "BasicDoubleton.h"
$LCOV_CMD --list coverage.info | grep "SharedDoubleton.h"
$LCOV_CMD --list coverage.info | grep "GenericJet.h"
$LCOV_CMD --list coverage.info | grep "DDEJetSection.h"

echo ">> Done."
