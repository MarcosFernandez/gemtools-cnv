# PROJECT: GEM-Tools-cnv library
# FILE: Makefile
# DATE: 01/05/2015
# AUTHOR(S): Marcos Fernandez <mfernandez@pcb.ub.es>
# DESCRIPTION: Builds the GT-tools-cnv
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

CFLAGS = -O2 -g -Wall
LDFLAGS = -lz -lm
SOURCES = mrcanavar.c utils.c prep.c sam.c callcnv.c gcnorm.c
OBJECTS = $(SOURCES:.c=.o)

GEM_TOOLS_CNV=mrcanavar
GEM_TOOLS_CNV_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS_CNV))


		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

all: GEM_TOOLS_CNV_FLAGS= $(CFLAGS) $(ARCH_FLAGS) $(LDFLAGS)
all: $(GEM_TOOLS_CNV_BIN)
	

static: GEM_TOOLS_CNV_FLAGS= $(CFLAGS) $(ARCH_FLAGS) $(LDFLAGS) -static
static: $(GEM_TOOLS_CNV_BIN)

debug: GEM_TOOLS_CNV_FLAGS= $(CFLAGS) $(ARCH_FLAGS) $(LDFLAGS) $(DEBUG_FLAGS)
debug: $(GEM_TOOLS_CNV_BIN)


$(GEM_TOOLS_CNV_BIN): $(GEM_TOOLS_CNV) 
	$(CC) -o $@ $(SOURCES) $(GEM_TOOLS_CNV_FLAGS)  



