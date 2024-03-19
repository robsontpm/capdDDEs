Nonrigorous, non-autonomous demo
================================

This example replicates the result of the 
[exemplary program]((http://capd.sourceforge.net/capdDynSys/docs/html/examples_RosslerChaoticDynamics.html) 
from the CAPD library about the existence of the attractor and the symbolic dynamic in the 
[Rossler system](http://www.scholarpedia.org/article/Rossler_attractor).

This example shows that the capdDDEs code can be used in a quite similar way to the ODE code from CAPD.
However, if you are going to work only with ODEs (not true DDEs), the CAPD is better for this task,
as it is a little easier, can integrate backward in time, and is much faster. 

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

    ./dde-vs-ode-code

- have a nice day!

NOTE
----

Th program ```cad-original-code.cpp``` is a sligthly cropped version
of the original example from CAPD library, just for reference.