# -*- Mode: Makefile -*-

# trace the chain of include makefiles
makefiles += NASA4_exec

# This is an embedded boundary example
USE_EB=TRUE

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = plasma

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../lib
##CHOMBO_HOME = /g/g17/massa5/Chombo-3.2/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTimeDependent EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTimeDependent AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ./src_CN \
           ./src_POTENTIAL \
           ./src_SPECIES \
           ./src_PHYSICS

# input file for 'run' target
INPUT = BL.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables



# application-specific targets
varsall:
	@echo "syslibflags=$(syslibflags)"
