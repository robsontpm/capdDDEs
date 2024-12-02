# this will make it into a separate subfolder of a given name, trailing / is important!) 
PACKAGE = demo-elninio/

# the defaut output folder (defined above) is ignored if those two are given (please keep trailing '/' !)
OUTDIR = ./bin/
OBJDIR = ./.obj/

# a list of all the programs to compile (from this directory), separated by single space.
# i.e. .cpp files that contain main() function.
#EXECUTABLES = nonrig timemap-rig poincare-rig 
EXECUTABLES = 1_nonrig 2_manipulations 3_rigorous_timemap 4_periodic_candidate

# a list of all your units to be linked with your programs (other .cpp files), separated by single space
# i.e. all other .cpp files that does not have main() inside.
OTHERS = setup

# 'run' the true Makefile
include ../makefile-run.mk