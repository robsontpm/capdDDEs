######################################################################
#          !!! this is version for bundled with CAPD !!!             #
######################################################################
#                                                                    #
# You need to set the paths to your compiled CAPD                    #
#                                                                    #
######################################################################

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
# Please uncomment one. 
CAPDSCRIPT := capd-config
#CAPDSCRIPT := mpcapd-config

# BUT REMEMBER: you need to update all makefile-template.mk files!

