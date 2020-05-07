LIBRARY = TH4D
OBJDIR = obj
DEPDIR = $(OBJDIR)/dep
SRCDIR = ./
INCDIR = ./
ADDITIONAL_ROOTMAPLIBRARY= -rml libHist.so

ifdef SKFlat_WD
  include $(SKFlat_WD)/makefile.common
else
  include makefile.common
endif

