# relevant targets in subdirectory makefiles
TOPTARGETS := all clean
SUBDIRS := $(wildcard */.)

# call 'make all' (same as 'make') or 'make clean' in all subdirectories
.PHONY: $(TOPTARGETS) $(SUBDIRS)
$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	make -C $@ $(MAKECMDGOALS)
