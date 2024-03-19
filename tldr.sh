#!/bin/bash

cd external

# add permission to helper scripts if  and run build

chmod a+x capd/configure
chmod a+x capd-build.sh
./capd-build.sh

# then, after capd is build without errors, then compile programs.
# First go back to root directory of this project

cd ..

# then make the selected program or demo eg:

make examples/mackey_glass_stable

# compiled programs will be in ./bin directory

cd bin/examples/mackey_glass_stable && ls

# sample output of the above ls commands:
#
#    find_and_prove_n6
#
# finally, run proofs or demos, e.g. :

./find_and_prove_n6

# you can see list of available programs and demos by running
# in the main directory the following command:

cd ../../..
make list

# if you are interested in more documentation, see ./docs