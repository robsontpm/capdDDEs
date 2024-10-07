#!/bin/bash

# this script build local version of CAPD
# we follow rougly http://capd.sourceforge.net/capdDynSys/docs/html/capd_compilation.html

# Working Directory is where this script is located
WD=`pwd`
BIN_DIR="$WD/../bin"
BUILD_DIR="$BIN_DIR/capd_build"
CAPD_SRC_DIR="$WD/capd"
# Path relative to BUILD_DIR, where configure script sits
CAPD_REL_PATH="../../external/capd"
# where to install the library (locally)
PREFIX="$BIN_DIR/capd"

# remove build to start anew
if [ -d "$BUILD_DIR" ]
then
	rm -rf "$BUILD_DIR"
fi

# make build and go there
mkdir -p "$BUILD_DIR"

if [ -d $PREFIX ]
then
	echo "$PREFIX exists"
else
	mkdir -p "$PREFIX"
fi
# those three lines are to get abs path
cd "$PREFIX"
PREFIX=`pwd`
cd "$WD"
echo "$PREFIX"

# make sure GCC and G++ are the compilers
# clang may make strange errors...
CXX=g++
CC=gcc

# do the compilation (outside CAPD directory)
echo "cd to build directory: $BUILD_DIR"
cd "$BUILD_DIR"
cmake --install-prefix "$PREFIX" "$CAPD_REL_PATH" -DCAPD_ENABLE_MULTIPRECISION=true
cmake --build .
cmake --install .

# copy the missing include file into installation directory
# bug in capd cmake files?
cp -r "$CAPD_SRC_DIR"/capdAlg/include/capd/intervals/intra "$PREFIX"/include/capd/intervals
