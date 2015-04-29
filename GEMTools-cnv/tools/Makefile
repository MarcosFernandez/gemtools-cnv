#==================================================================================================
# PROJECT: GEM-Tools-cnv library
# FILE: Makefile
# DATE: 01/05/2015
# AUTHOR(S): Marcos Fernandez <mfernandez@pcb.ub.es>
# DESCRIPTION: Builds the GT-tools-cnv
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GEM_TOOLS_CNV=kmermaker

GEM_TOOLS_CNV_SRC=$(addsuffix .c, $(GEM_TOOLS_CNV))
GEM_TOOLS_CNV_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS_CNV))
LIBS:=-lm


ifeq ($(HAVE_ZLIB),1)
LIBS:=$(LIBS) -lz
endif
ifeq ($(HAVE_BZLIB),1)
LIBS:=$(LIBS) -lbz2
endif

all: GEM_TOOLS_CNV_FLAGS=-O4 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(GEM_TOOLS_CNV_BIN)
	
static: GEM_TOOLS_CNV_FLAGS=-O4 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -static
static: $(GEM_TOOLS_CNV_BIN)

debug: GEM_TOOLS_CNV_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(GEM_TOOLS_CNV_BIN)

$(GEM_TOOLS_CNV_BIN): $(GEM_TOOLS_CNV_SRC)
	$(CC) $(GEM_TOOLS_CNV_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(EXTRA_LIBS)
	
                            
