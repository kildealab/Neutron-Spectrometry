######################################################################
#
#	Makefile for Georges_test 
#
#	J. Kildea  28 July, 2016 
#
#######################################################################

#########################################################################
# ROOT stuff
#########################################################################
ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs)
ROOTGLIBS    = $(shell root-config --glibs)

LIBS            = $(ROOTLIBS)
GLIBS           = $(ROOTGLIBS)

OAWGDEFFLAGS = -D_OAWG -D_OAWG_DIR="\"`pwd`\"" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS

#########################################################################
# General stuff
#########################################################################
CPP = g++
CFLAGS = -o $@ -Wall -O -g $(ROOTCFLAGS) 

ALLLIBS = $(LIBS) $(GLIBS) 

#########################################################################
# Targets
#########################################################################
all: main.exe

main.exe: main.o
	$(CPP) $(CFLAGS) main.o $(ALLLIBS)

main.o: main.cpp 
	$(CPP) -c $(CFLAGS) main.cpp

# tidy up
clean: 
	rm -rf *.o main.exe
	
