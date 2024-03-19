Utility programs
================

This folder contains programs that might be of 
use to a wider audience. The programs are as follows:

 - benchmark

   A simple program that shows difference between
   non-expanding and expanding repreentations
   in terms of the accuracy. See FoCM 2023 (2024) paper
   for more details. 

 - converter

   a set of programs to convert between various
   file types used in CAPD and CAPDDDEs packages.
   See documentation inside the folder for more details. 

 - tests

   contains a set of tests, that should provide
   strong evidence the code is working as expected. 
   **If the tests fails** on your system, please
   report this to the developer, providing as much
   information as you can, including your computer
   specyfications (especially model of the processor).

To compile programs, invoke ```make utils/program-name```
e.g. ```make utils/converter```. The executables go
directly to **bin** folder in the root directory of 
the repository. 