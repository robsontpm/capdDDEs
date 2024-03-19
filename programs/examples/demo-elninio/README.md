Nonrigorous, non-autonomous demo
================================

The equation used in this example comes from works:

[^1]: Delayed Feedback Versus Seasonal Forcing: Resonance Phenomena in an El Nin͂o Southern Oscillation Model
  Andrew Keane, Bernd Krauskopf, and Claire Postlethwaite, (2015), https://doi.org/10.1137/140998676
  
[^2]: The effect of state dependence in a delay differential equation model for the El Niño Southern Oscillation
  Andrew Keane, Bernd Krauskopf and Henk A. Dijkstra, (2019), https://doi.org/10.1098/rsta.2018.0121

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

    ./draw-nonrig

- to see possible options of the program, type

    ./draw-nonrig --help

- you may also try to run 'run-demos.sh' from this demo root directory:

    bash run-demos.sh

- have a nice day!