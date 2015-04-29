CPP      = g++
CPPFLAGS   = -Wall -O3 -march=native
CPPFLAGS_C  = -Wall -O2


all: galculator

simple: galculator-comp

galculator:	galculator.cpp
	$(CPP) $(CPPFLAGS) -o ../bin/galculator galculator.cpp

galculator-comp:	galculator.cpp
	$(CPP) $(CPPFLAGS_C) -o  ../bin/galculator galculator.cpp
