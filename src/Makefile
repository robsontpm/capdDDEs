include ../build-conf.mk

include ../build-common.mk

VPATH = $(shell head -n 1 libvpath.txt)

OTHERS = $(shell head -n 1 libcapdddes.txt)

# directory where object and dependancy files will be created
OBJDIR = ${PROJECT_PATH}/.obj/
OUTDIR = ${PROJECT_PATH}/bin/

# other setting important for this and for the above are dfined 
# in the main build-conf.mk file (see root directory)

OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
OBJ_FILES = ${OTHERS_OBJ} 

default: 
	@echo "This is config test"
	@echo "==================="
	# setting compiler and linker flags
	@echo "CXX:        ${CXX}"
	@echo "OTHERLIBS:  ${OTHERLIBS}"
	@echo "CAPDFLAGS:  ${CAPDFLAGS}"
	@echo "CAPDLIBS:   ${CAPDLIBS}"
	@echo "CXXFLAGS:   ${CXXFLAGS}"
	@echo "OBJDIR:     ${OBJDIR}"
	@echo "OUTDIR:     ${OUTDIR}"
	@echo "OTHERS:     ${OTHERS}"
	@echo "OTHERS_OBJ: ${OTHERS_OBJ}"
	@echo "OBJ_FILES:  ${OBJ_FILES}"
	@echo "VPATH:      ${VPATH}"
	@echo "==================="
	@echo "call make all to build a program"

all:  ${OTHERS_OBJ}

## rule to link executables
#${PROGS}: % : ${OBJDIR}%.o ${OTHERS_OBJ}
#	@mkdir -p ${OUTDIR}
#	${CXX} ${LINKERFLAGS} -o ${OUTDIR}$@ $< ${OTHERS_OBJ} ${CAPDLIBS}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}

#rule to compile .cpp files and generate corresponding files with dependencies
${OBJ_FILES}: ${OBJDIR}%.o : %.cpp
	@mkdir -p ${OBJDIR}
	$(CXX) ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<

# rule to clean all object files, dependencies and executables
.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o ${OBJDIR}*.o.d





