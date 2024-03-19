###################################################################
# this is a facade that redirects to src/programs/[given_program] #
###################################################################

include ./build-conf.mk

# Tell make that these are phony targets
.PHONY: all list clean default test

# TODO: (important!) make sure demos have different target-name than examples
# TODOL (important!) now I need to make sure the names do not collide...

# list all programs in directory (only names, not the whole paths)
RESULTS=$(shell find ./programs/results/ -maxdepth 1 -mindepth 1 -type d -printf "results/%P\n" )
EXAMPLES=$(shell find ./programs/examples/ -maxdepth 1 -mindepth 1 -type d -printf "examples/%P\n" )
DEVEL=$(shell find ./programs/devel/ -maxdepth 1 -mindepth 1 -type d -printf "devel/%P\n" )
UTILS=$(shell find ./programs/utils/ -maxdepth 1 -mindepth 1 -type d -printf "utils/%P\n" )
export ${RESULTS}
export ${EXAMPLES}
export ${DEVEL}
export ${UTILS}
PROGRAMS=${RESULTS} ${EXAMPLES} ${DEVEL} ${UTILS}
export ${PROGRAMS}

#@ at the begining tells make not to print the command
list:	
	@echo "=== List of possible make targets ==="
	@echo "(type make <target/name>) to compile"
	@echo "1) examples (simple progs for learning):" 
	@for item in ${EXAMPLES}; do echo "    * $$item"; done
	@echo "2) results (from publications):" 
	@for item in ${RESULTS}; do echo "    * $$item"; done
	@echo "3) devel (in development):" 
	@for item in ${DEVEL}; do echo "    * $$item"; done	
	@echo "4) utils (for general use):" 
	@for item in ${UTILS}; do echo "    * $$item"; done		
	@echo "5) misc:" 
	@echo "    * all\n    * test\n    * clean\n    * list\n    * programs"
	@echo "==================================="

$(PROGRAMS):
	@echo "Building program: $@..."
	make -C ./programs/$@ all

default: list

# Clean up the executable files
clean:
	rm -rf ./.obj
	
programs: $(PROGRAMS)

test:
	make -C utils/tests all

all: programs

#lib: capdddes.a
#
#capdddes.a: 

