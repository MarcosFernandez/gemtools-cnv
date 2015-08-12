# PROJECT: GEM-Tools-cnv library
# FILE: Makefile
# DATE: 01/05/2015
# AUTHOR(S): Marcos Fernandez <mfernandez@pcb.ub.es>
# DESCRIPTION: Builds the GT-tools-cnv
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

CFLAGS = -c -O2 -g -Wall
LDFLAGS = -lz -lm
SOURCES = mrcanavar.c utils.c prep.c sam.c callcnv.c gcnorm.c
OBJECTS = $(SOURCES:.c=.o)

EXECUTABLE = mrcanavar
BIN = $(addprefix $(FOLDER_BIN)/, $(EXECUTABLE))

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
	cp mrcanavar $(BIN)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ \#* 





