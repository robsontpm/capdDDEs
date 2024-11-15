######################################################################
#        A Common file to setup common symbols used in Makefiles     #
######################################################################
#               NOTICE: You not need to modify this file.            #
######################################################################

# ?= is used to allow user to evrload those in their Makefile
GPROF_OPTION ?=
GGDB_OPTION ?=

# setting compiler and linker flags
#OTHERLIBS = -lboost_iostreams -lboost_system -lboost_filesystem ${GPROF_OPTION}
OTHERLIBS ?= ${GPROF_OPTION}
CAPDFLAGS = `${CAPDBINDIR}${CAPDSCRIPT}  --cflags`
CAPDLIBS = `${CAPDBINDIR}${CAPDSCRIPT}  --libs` ${OTHERLIBS}
CXXWARNINGFLAGS ?=
CXXEXTRAFLAGS ?= -O2 -std=c++17
LINKERFLAGS ?= ${EXTRALINKEROPTIONS} -L${CAPDLIBDIR}
CXXFLAGS += ${CAPDFLAGS} ${CXXEXTRAFLAGS} ${CXXWARNINGFLAGS} -I${PROJECT_PATH}/include ${GPROF_OPTION} ${GGDB_OPTION}


