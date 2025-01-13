A full demo
===========

This demo was prepared during a stay in SZTE Bolyai Instutute, Hungary.

The examples should be inspected in the order: ``1_nonrig.cpp``
then ``2_manipulations``, etc.,
as the later programs does not contain comments explained already in the previous files. 

The programs start with a numeral to keep the order of files in the operating systems.

Extra files generated from programs will have names that will put them behind any other files. 
Please do not run the programs from the source directory (e.g. invoking ```./bin/1_nonrig```),
because then the files will appear in the source code directory and it will be a mess. 

TODO: Expand the documentation.

Compilation
-----------

- before compiling this demo, you need to compile extrenal libraries
  first, go to the root directory fo the project and read [README.md](../../../README.md)
  there.

- When you have compiled the external libraries, you can return here 
  and compile this demo 

- to compile this demo you need to invoke in this directory:

    make all
    
- it will create ./bin subdirectory

    cd ./bin

- to run the program, simply type e.g. :

    ./1_nonrig
    