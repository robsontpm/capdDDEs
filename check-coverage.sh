#!/bin/bash

WD=$(pwd)
BUILD_DIR="build"

# Configure with Coverage enabled
cd $BUILD_DIR
echo ">> Configuring CMake with ENABLE_COVERAGE=ON..."
cmake -DENABLE_COVERAGE=ON -DBUILD_PROGRAMS=OFF -DCAPD_DIR="$CAPD_DIR" ..

# make all tests
make -j4

echo ">> Running ALL tests..."
ctest -V -R

# Generate Coverage Report
echo ">> Generating Coverage Report..."
# Capture coverage data using lcov
# Ignoring mismatch errors which can happen with templates/headers
lcov --capture --directory . --output-file coverage.info --ignore-errors gcov --base-directory ..

# Generate summary for all current tests
echo ">> Coverage for specyfic components:"
lcov --list coverage.info | grep "DoubletonInterface.h"
lcov --list coverage.info | grep "BasicDoubleton.h"
lcov --list coverage.info | grep "SharedDoubleton.h"
lcov --list coverage.info | grep "GenericJet.h"

echo ">> Done."
