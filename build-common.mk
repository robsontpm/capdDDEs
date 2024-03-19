######################################################################
#        A Common file to setup common symbols used in Makefiles     #
######################################################################
#               NOTICE: You not need to modify this file.            #
######################################################################

# I had to override it here, as Eclipse builder do not accept LLMV 
# toolchain, and override it with GCC...
# TODO: my current clang does not compile with mpcapd (Multiprec) 
#export CXX=clang++
#export CC=clang
export CXX=g++
export CC=gcc

# choose one of those 
#   first: support ggdb, slower execution time; 
#   second: no debuger info, faster execution;
#GGDB_OPTION=-ggdb3 -g
GGDB_OPTION=

# choose one of those
#   first:  gprof program (profiling, speed testing), slower execution time; 
#   second: no profiler info, faster execution;
#GPROF_OPTION=-pg
GPROF_OPTION=

# setting compiler and linker flags
#OTHERLIBS = -lboost_iostreams -lboost_system -lboost_filesystem ${GPROF_OPTION}
OTHERLIBS = ${GPROF_OPTION}
CAPDFLAGS = `${CAPDBINDIR}${CAPDSCRIPT}  --cflags`
CAPDLIBS = `${CAPDBINDIR}${CAPDSCRIPT}  --libs` ${OTHERLIBS}
#CXXWARNINGFLAGS = -pedantic -Wall -Wclobbered -Wempty-body -Wignored-qualifiers -Wmissing-field-initializers -Wsign-compare -Wtype-limits -Wuninitialized -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual
#CXXWARNINGFLAGS = -pedantic -Wall -Wclobbered -Wempty-body -Wignored-qualifiers -Wmissing-field-initializers -Wsign-compare -Wtype-limits -Wuninitialized -Wundef -Wcast-align -Wwrite-strings -Wlogical-op
CXXWARNINGFLAGS = 
CXXEXTRAFLAGS = -O2 -std=c++17
# static is used to prevent linker error for CAPD. Other methods to avoid it are described in CAPD docs. TODO: (NOT URGENT): add better linking (dynamic?) 
LINKERFLAGS = -static -L${CAPDLIBDIR}
#LINKERFLAGS = 
CXXFLAGS += ${CAPDFLAGS} ${CXXEXTRAFLAGS} ${CXXWARNINGFLAGS} -I${PROJECT_PATH}/include ${GPROF_OPTION} ${GGDB_OPTION}


