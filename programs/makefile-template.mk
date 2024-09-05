#######################################################################
#                                                                     #
# This is a template for any Makefile in 'programs                    #
#                                                                     #
#######################################################################
#                                                                     #
# concrete templates would be made for each of the subfolder          #
# see the documentation in those files to learn how to configure      #
# your own Makefiles for programs                                     #
#                                                                     #
#######################################################################

# directory where object and dependancy files will be created
# ?= is used to allow user to evrload those in their Makefile
OBJDIR ?= ${PROJECT_PATH}/.obj/${PROGRAMS_SUBDIR}${PACKAGE}
OUTDIR ?= ${PROJECT_PATH}/bin/${PROGRAMS_SUBDIR}${PACKAGE}

# other setting important for this and for the above are defined 
# in the main build-conf.mk file (see root directory)

# other object .o files from the capdddes lib
# (added here, as I sometimes develop both lib and the programs to 
#  test them at the same time, and to not forget that i need to 
#  recompile them)
LIBITEMS = $(shell head -n 1 ${PROJECT_PATH}/src/libcapdddes.txt)
VPATHITEMS = $(shell head -n 1 ${PROJECT_PATH}/src/libvpath.txt)
VPATH = ${VPATHITEMS:%=${PROJECT_PATH}/src/%}

# now we list all OTHERS objects (user defined), EXECUTABLES 
# (user defined) and LIBS objects (compiled items from the ddes lib)
# see note above why I do it tat way
OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
LIBS_OBJ = ${LIBITEMS:%=${PROJECT_PATH}/.obj/%.o}

# put all PROGRAM object files (with paths) to this variable
# we will treat LIBS_OBJ separately
OBJ_FILES = ${OTHERS_OBJ} ${EXECUTABLES:%=${OBJDIR}%.o}

# running just make will print this info
default: 
	@echo "This is config test"
	@echo "==================="
	# setting compiler and linker flags
	@echo "CXX:         ${CXX}"
	@echo "OTHERLIBS:   ${OTHERLIBS}"
	@echo "CAPDFLAGS:   ${CAPDFLAGS}"
	@echo "CAPDLIBS:    ${CAPDLIBS}"
	@echo "CXXFLAGS:    ${CXXFLAGS}"
	@echo "OBJDIR:      ${OBJDIR}"
	@echo "OUTDIR:      ${OUTDIR}"
	@echo "EXECUTABLES: ${EXECUTABLES}"
	@echo "OTHERS:      ${OTHERS}"
	@echo "OTHERS_OBJ:  ${OTHERS_OBJ}"
	@echo "OBJ_FILES:   ${OBJ_FILES}"
	@echo "VPATH:       ${VPATH}"
	@echo "==================="
	@echo "call 'make all' to build the program(s)"
	@echo "call 'make rmexecutables' to remove program(s)"
	@echo "call 'make clean' to remove build .o/.d files"	

# running 'make all' will try to compile all EXECUTABLES
# but they are dependent on OTHERS and LIBS_OBJ 
all: ${EXECUTABLES}

# rule to link executables (with OTHERS and also LIBS from capdddes lib)
# NOTE: $@ - means full target name, i.e. program-name, 
# NOTE: $< - means first dependency, i.e. program-name.o
${EXECUTABLES}: % : ${OBJDIR}%.o ${LIBS_OBJ} ${OTHERS_OBJ}
	@mkdir -p ${OUTDIR}
	${CXX} ${LINKERFLAGS} -o ${OUTDIR}$@ $< ${LIBS_OBJ} ${OTHERS_OBJ} ${CAPDLIBS}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}
-include ${LIBS_OBJ:%=%.d}

#rule to compile .cpp files and generate corresponding files with dependencies
${OBJ_FILES}: ${OBJDIR}%.o : %.cpp
	@mkdir -p ${OBJDIR}
	$(CXX) ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<
	
# TODO: this is brzydkie, trzebaby to zamienić	
${LIBS_OBJ}: ${PROJECT_PATH}/.obj/%.o : %.cpp
	@mkdir -p ${PROJECT_PATH}/.obj/
	$(CXX) ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<

# rule to clean all object files, dependencies and executables
.PHONY: clean rmexecutables ${EXTRA_PHONY}
clean:
	rm -f ${OBJDIR}*.o ${OBJDIR}*.o.d
	
rmexecutables:
	echo "${EXECUTABLES}" | xargs -d ' ' -I{} rm -f ${OUTDIR}{}

