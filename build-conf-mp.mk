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
CAPDBINDIR := $(PROJECT_PATH)/bin/capd/bin/
CAPDLIBDIR := $(PROJECT_PATH)/bin/capd/lib/

# this is to tell to use multiprecision-enabled CAPD
CAPDSCRIPT := mpcapd-config
