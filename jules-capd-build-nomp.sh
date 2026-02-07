#!/bin/bash
# Adjusted for running from repo root without MPFR

WD=$(pwd)
BIN_DIR="$WD/bin"
BUILD_DIR="$BIN_DIR/capd_build"
CAPD_SRC_DIR="$WD/external/capd"
CAPD_REL_PATH="$WD/external/capd"

if [ -d "$BUILD_DIR" ]; then rm -rf "$BUILD_DIR"; fi
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake "$CAPD_REL_PATH" -DCAPD_ENABLE_MULTIPRECISION=OFF
make -j4
