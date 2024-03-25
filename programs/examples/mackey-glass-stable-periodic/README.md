Rigorous computer assisted proof of a periodic orbit
====================================================

This example replicates the result of the first FoCM paper 
(see references in [README.md](../../../README.md)),
but it is written in a more concise way. It shows usage of
the helper classes.

Compilation
-----------

 * before compiling this demo, you need to compile extrenal libraries
   first, go to the root directory fo the project and read [README.md](../../../README.md)
   there.
 
 * When you have compiled the external libraries, you can return here 
   and compile this demo 
 
 * to compile this demo you need to invoke in this directory:
 
     make all
     
 * it will create ``./bin`` subdirectory
 
     ls ./bin
 
 * the computation is split into 3 steps:
    1. Finding the candidate nonrigorously, done by ``nonrig-find.cpp`` 
       It finds also the apparently good coordinates for the initial set on the section.
       The set will have form ``X = x + C * r0``, where ``x`` is the good approfimation of the ``p_0``
       segment of the periodic solution ``p : R \to R``. We also compute good section on which ``X`` is given.
    2. Rigorously computing the inverse of the matrix ``C``, so that the ``P(X)`` can be later 
       transformed into good coordinates. This step is done with the help of general set of 
       programs available in ``programs/utils/converter``    
    3. Computing ``P(X)`` and comparing ``C^-1 (P(X) - x0) =: Pr0`` to ``r0``.
	   Here, ``P`` is the Poincare map that comes back to the section on which ``X`` was defined. 
	   This part is done by ``rig-prove.cpp``.
       
 * all the setps can be run in order, simply by running:
 
	bash run-proof.sh
    
 * We are using Schauder Fixed Point theorem to infer that if ``Pr0 \subset r0``, then there exists
   ``p \in X`` that is a periodic point of ``P``, thus a periodic solution to the DDE. 
   For otherassumptions of Schauder Fixed Point theorem, the compactness comes from the smoothing
   of solutions in DDEs, and the set ``X`` is convex by definition. See papers in references for more details 
   [README.md](../../../README.md)
 
 * have a nice day!

NOTE
----

You can try to change the parameters, but the program will not work
on apparently unstable periodic solutions. It might also
fail on a more complicated solutions, e.g. period-2 points for 
Poincare map.

How to handle more complicated dynamics, see our other papers
mentioned in main [README.md](../../../README.md) and references 
therein.