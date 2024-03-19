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

# for some reason I need to add permissions to execute to many files in this dir
chmod -R a+x "$CAPD_SRC_DIR/capdMake/libcapd/"

# make sure GCC and G++ are the compilers
# clang may make strange errors...
CXX=g++
CC=gcc

# do the compilation (outside CAPD directory)
echo "cd to build directory: $BUILD_DIR"
cd "$BUILD_DIR"
chmod a+x "$CAPD_REL_PATH/configure"
$CAPD_REL_PATH/configure --prefix "$PREFIX" --with-mpfr=yes
# FOR SOME REASON I need to use relative path in the configure invocation
# otherwise I get strange errors from the CAPD configuration process...

# build it (it should go to $PREFIX)
make -j 8

# move everything to bin/capd for future use (dir.conf uses this location)
make install
