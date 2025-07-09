Command Line Parser - Helper class
==================================

This is an exemplary program that uses the simple 
command line parser shipped with the library. 

I use this parser a lot in my programs, to be able to
pass variables in a convenient way, without relying
on modern and very complicated external libraries. 

The command line parser also has the ability to 
work with capd::intervals, so e.g. it parses
non-representable floats into small-diameter intervals. 

It also has the ability to output the information 
on the computation, i.e. the command line used
and the date and time when the program was run. 
I usually put those in the headers (on in the end)
of output files to guide me in my reasearch after 
returning to the experiments after few weeks/months. 
This is a great help! Trust me!

NOTE
----

I am thinking about moving this parser outside,
to a completely separate repository. This is 
a TODO: for future.

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

- to run the program, simply type

    ./cmd-parser

- to see possible options of the program, type

    ./cmd-parser --help

- have a nice day!