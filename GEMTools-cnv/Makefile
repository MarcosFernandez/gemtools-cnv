#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/05/2015
# AUTHOR(S): Marcos Fernandez Callejo <mfernandez@pcb.ub.es>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)
include Makefile.mk

all: release

static:	setup
#	$(MAKE) --directory=src
#	$(MAKE) --directory=resources/src
	$(MAKE) --directory=tools static
	$(MAKE) --directory=mrcanavar static
	$(MAKE) --directory=bam2fastq/samtools static
	$(MAKE) --directory=bam2fastq static
	$(MAKE) --directory=galculator static
	$(MAKE) --directory=src static


release: setup	
#	$(MAKE) --directory=src
#	$(MAKE) --directory=resources/src
	$(MAKE) --directory=tools
	$(MAKE) --directory=mrcanavar
	$(MAKE) --directory=bam2fastq/samtools
	$(MAKE) --directory=bam2fastq
	$(MAKE) --directory=galculator
	$(MAKE) --directory=src

debug: setup
#	$(MAKE) --directory=src   debug
#	$(MAKE) --directory=resources/src   debug
	$(MAKE) --directory=tools debug
	$(MAKE) --directory=mrcanavar debug
	$(MAKE) --directory=bam2fastq/samtools
	$(MAKE) --directory=bam2fastq debug
	$(MAKE) --directory=galculator debug
	$(MAKE) --directory=src debug

check: setup debug
	$(MAKE) --directory=test  check
	
setup: 
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)

clean:
#	$(MAKE) --directory=test  clean
	@rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	$(MAKE) --directory=bam2fastq clean
        $(MAKE) --directory=bam2fastq/samtools  clean
