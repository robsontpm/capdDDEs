#!/bin/bash

# download CAPD external library first
git submodule update --init

# go to the folder where CAPD library is
cd external

# this is no longer needed because of a git submodule
# instead we do submodule --init above.
# chmod a+x capd/configure

# add permission to helper scripts if  and run build
chmod a+x capd-build.sh
./capd-build.sh

# then, after capd is build without errors, then compile programs.
# First go back to root directory of this project

cd ..

# then make the selected program or demo eg:

make examples/rossler-ode-vs-dde-code

# compiled programs will be in ./bin directory

cd programs/examples/rossler-ode-vs-dde-code/bin && ls

# sample output of the above ls commands:
#
#    capd-original-code  dde-vs-ode-code  nonrig-plotter
#
# finally, run proofs or demos, e.g. :

./dde-vs-ode-code

# you can see list of available programs and demos by running
# in the main directory the following command:

cd ../../../..
make list

# if you are interested in more documentation, see ./docs