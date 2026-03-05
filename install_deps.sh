#!/bin/bash
set -e
mkdir -p ${HOME}/deps
cd ${HOME}/deps

# Install Boost 1.84.0 from sourceforge
if [ ! -d "${HOME}/deps/include/boost" ]; then
    rm -f boost_1_84_0.tar.gz*
    wget -O boost_1_84_0.tar.gz "https://sourceforge.net/projects/boost/files/boost/1.84.0/boost_1_84_0.tar.gz/download"
    tar -xzf boost_1_84_0.tar.gz
    cd boost_1_84_0
    ./bootstrap.sh --with-libraries=test
    ./b2 install --prefix=${HOME}/deps
    cd ..
fi
