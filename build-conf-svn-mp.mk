######################################################################
#          !!! this is version for your own CAPD !!!                 #
######################################################################
#                                                                    #
# You need to set the paths to your compiled CAPD                    #
#                                                                    #
######################################################################

# I had to override it here, as Eclipse builder do not accept LLMV 
# toolchain, and override it with GCC...
# TODO: my current clang does not compile with mpcapd (Multiprec) 
#export CXX=clang++
#export CC=clang
export CXX=g++
export CC=gcc

# those are taken from some stackoverflow and will return
# the paths to __this__ file (..._FILEPATH) and the absolute path
# to the folder (..._DIRPATH). It will work even for Makefiles
# where I include the build-conf.mk file. 
MAKEFILE_ABS_FILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MAKEFILE_ABS_DIRPATH := $(patsubst %/,%,$(dir $(MAKEFILE_ABS_FILEPATH)))
PROJECT_PATH := $(MAKEFILE_ABS_DIRPATH)

# if no CAPDBINDIR and CAPDLIBDIR is given, then a shell global
# variable will be used. Therefore, if you setup your shell (bash)
# variables to your external CAPD compilation,  you can use it
# to compile 
# 
# just run something like that in shell: 
#
#    export CAPDBINDIR=~/eclipse-workspace/repo/capd_svn_install/bin/
#    export CAPDLIBDIR=~/eclipse-workspace/repo/capd_svn_install/lib/
#
# or add such lines into your .bashrc or similar profile file

# this is to tell which version of CAPD to use
# (either capd, or mpcapd with multiprecision)
CAPDSCRIPT := mpcapd-config

# static is used to prevent linker error for CAPD. 
# Other methods to avoid it are described in CAPD docs. 
#EXTRALINKEROPTIONS=-static
# uncomment for: support ggdb, slower execution time; 
#GGDB_OPTION=-ggdb3 -g
# uncomment for: gprof program (profiling, speed testing), slower execution time; 
#GPROF_OPTION=-pg
# TODO: think if uncomment some of those to get cleaner code?
#CXXWARNINGFLAGS = -pedantic -Wall -Wclobbered -Wempty-body -Wignored-qualifiers -Wmissing-field-initializers -Wsign-compare -Wtype-limits -Wuninitialized -Wundef -Wcast-align -Wwrite-strings -Wlogical-op
