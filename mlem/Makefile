######################################################################
#
#	Makefile for neutron unfolding program: main.cpp
#
#######################################################################

#########################################################################
# ROOT libraries & flags
#########################################################################
ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs)
ROOTGLIBS    = $(shell root-config --glibs)

LIBS            = $(ROOTLIBS)
GLIBS           = $(ROOTGLIBS)

OAWGDEFFLAGS = -D_OAWG -D_OAWG_DIR="\"`pwd`\"" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS

#########################################################################
# General libraries & flags
#########################################################################
CPP = g++

GIT_COMMIT := $(shell git describe --abbrev=7 --dirty --always --tags)
GIT_CFLAG = -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CFLAGS = -o $@ -Wall -O -g -std=c++11 $(GIT_CFLAG)

LFLAGS = -Wall -O -g $(ROOTCFLAGS) 

OBJS = main.o fileio.o handle_args.o root_helpers.o physics_calculations.o

ALLLIBS = $(LIBS) $(GLIBS) 

#########################################################################
# Targets
#########################################################################
all: main.exe

main.exe: $(OBJS)
	$(CPP) $(LFLAGS) $(OBJS) $(ALLLIBS) -o main.exe

main.o: main.cpp fileio.h handle_args.h root_helpers.h physics_calculations.h
	$(CPP) -c $(CFLAGS) main.cpp

fileio.o: fileio.h fileio.cpp
	$(CPP) -c $(CFLAGS) fileio.cpp

handle_args.o: handle_args.h handle_args.cpp
	$(CPP) -c $(CFLAGS) handle_args.cpp

root_helpers.o: root_helpers.h root_helpers.cpp
	$(CPP) -c $(CFLAGS) $(ROOTCFLAGS) root_helpers.cpp

physics_calculations.o: physics_calculations.h physics_calculations.cpp
	$(CPP) -c $(CFLAGS) physics_calculations.cpp

# tidy up
clean: 
	rm -rf *.o main.exe
	
