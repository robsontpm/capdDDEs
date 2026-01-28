#!/bin/bash
set -e

# Define directories
BUILD_DIR="build-coverage"
REPORT_DIR="coverage_report"

# 1. Environment Setup
echo "Creating build directory..."
if [ -d "$BUILD_DIR" ]; then
    rm -rf "$BUILD_DIR"
fi
mkdir -p "$BUILD_DIR"

# 2. Configuration & Compilation
echo "Configuring and compiling..."
cd "$BUILD_DIR"
cmake -DENABLE_COVERAGE=ON ..
make -j$(nproc)

# 3. Counter Reset
echo "Zeroing counters..."
lcov --zerocounters --directory .

# 4. Test Execution
echo "Running tests..."
# We allow tests to fail but still want to generate the report
ctest --output-on-failure || echo "Some tests failed, proceeding to coverage report..."

# 5. Data Capture
echo "Capturing coverage..."
lcov --capture --directory . --output-file coverage.info --ignore-errors mismatch,gcov

# 6. Data Filtering
echo "Filtering coverage..."
# Remove system headers, tests, and external libraries (boost is covered by /usr/*)
lcov --remove coverage.info '/usr/*' '*/tests/*' '*/external/*' --output-file coverage_filtered.info --ignore-errors mismatch,gcov,unused

# 7. Report Generation
echo "Generating HTML report..."
genhtml coverage_filtered.info --output-directory "../$REPORT_DIR"

echo "Coverage report generated in $REPORT_DIR/index.html"
