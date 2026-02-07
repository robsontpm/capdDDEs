#!/bin/bash

# Adjusted for running from repo root

WD=$(pwd)
BIN_DIR="$WD/bin" # Build inside repo/bin
BUILD_DIR="$BIN_DIR/capd_build"
CAPD_SRC_DIR="$WD/external/capd"
CAPD_REL_PATH="$WD/external/capd" # Use absolute path to avoid confusion
PREFIX="$BUILD_DIR"

if [ -d "$BUILD_DIR" ]; then
	rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"

if [ ! -d $PREFIX ]; then
	mkdir -p "$PREFIX"
fi

cd "$PREFIX"
PREFIX=$(pwd)
cd "$WD"
echo "PREFIX: $PREFIX"

CXX=g++
CC=gcc

echo "cd to build directory: $BUILD_DIR"
cd "$BUILD_DIR"

cmake "$CAPD_REL_PATH" -DCAPD_ENABLE_MULTIPRECISION=true
make -j4
