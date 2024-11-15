######################################################################
#          !!! this is version for bundled with CAPD !!!             #
######################################################################
#                                                                    #
# this version will work if the capd is compiled from within         #
# './external/' directory in this project. Please check out if this  #
# directory contains 'capd/' folder, if it is so, then               #
# please build CAPD first (see README in 'external/' directory)      #
#                                                                    #
# NOTICE: You not need to modify this file.                          #
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

# this CAPD_BIN_DIR will be used in Makefiles for 
# programs and other executables
# directory under $CAPDBINDIR contains capd-config and mpcapd-config 
# scripts that are used in custom Makefiles to generate
# neccessary C++ compile/link flags.
# see http://capd.ii.uj.edu.pl for more info 
CAPDBINDIR := $(PROJECT_PATH)/bin/capd_build/bin/
CAPDLIBDIR := $(PROJECT_PATH)/bin/capd_build/lib/

# this is to tell which version of CAPD to use
# we use standard
CAPDSCRIPT := capd-config

# static is used to prevent linker error for CAPD. 
# Other methods to avoid it are described in CAPD docs. 
EXTRALINKEROPTIONS=-static
# uncomment for: support ggdb, slower execution time; 
#GGDB_OPTION=-ggdb3 -g
# uncomment for: gprof program (profiling, speed testing), slower execution time; 
#GPROF_OPTION=-pg
# TODO: think if uncomment some of those to get cleaner code?
#CXXWARNINGFLAGS = -pedantic -Wall -Wclobbered -Wempty-body -Wignored-qualifiers -Wmissing-field-initializers -Wsign-compare -Wtype-limits -Wuninitialized -Wundef -Wcast-align -Wwrite-strings -Wlogical-op
