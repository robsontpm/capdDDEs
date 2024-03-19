#######################################################################
#                                                                     #
# This is a template for Makefile in 'programs/devel/' directory      #
#                                                                     #
#######################################################################
#                                                                     #
# It should be used in the following way:                             #
#                                                                     #
#  * create program directory as ./[program_name]                     #
#  * create Makefile in the program directory                         #
#  * put the following in it:                                         #
#                                                                     #
#      EXECUTABLES = [program_name] (without cpp, may be more names)  #
#      OTHERS = names of other compile items (without cpp)            #
#      PACKAGE = [program_name]/                                      #
#      include ../makefile-run.mk                                     #
#                                                                     #
#  * now you can compile the program with the Makefile inside its     #
#    directory under ./[program_name] or directly from                #
#    Makefile in the project's main directory. Just call:             #
#                                                                     #
#        make [program_name]                                          #
#                                                                     #
#    e.g.                                                             #
#                                                                     #
#        make demo-elninio                                            # 
#                                                                     #
# Note 1: in [program_name] parentheses are part of the meta-name     #
# Note 2: PACKAGE must end with '/' and no trailing space !!!         #
#                                                                     # 
#######################################################################

# this will configure the makefile-template from the parent folder
PROGRAMS_SUBDIR=devel/

# Warning: the paths below are relative to the subfolder in this file's
#           folder as they will be invoked (inlcuded) from there!

# this will include helper file with paths configurable in one place
# and different configurations for CAPD
include ../../../build-conf.mk

# this will use symbols defined in the above file
include ../../../build-common.mk

# this will overwrite symbols from above file
# to configure it for needs of this specyfic kind of examples.
include ../../makefile-template.mk


