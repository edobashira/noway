# 
# Makefile for noway
# This needs the GNU version of make 
#

# When updating VERSION here, also update
# Version, Subversion and Patchlevel
# in NW_scheduler.h.
VERSION = 2.9.4

OSTYPE  = $(shell ./OSType)
BINTYPE	= $(shell ./BINType)
OBJDIR	= obj.$(BINTYPE)
DEPDIR  = .DEPEND.$(BINTYPE)
TDBDIR	= .TEMPLATE.$(BINTYPE)
$(shell mkdir -p $(OBJDIR))
$(shell mkdir -p $(DEPDIR))
$(shell mkdir -p $(TDBDIR)/Templates.DB)

prefix	= /share/spandh/packages/sprachworks
exec-prefix = $(prefix)/$(BINTYPE)
BINDIR	= $(exec-prefix)/bin

# g++
 CXX	= g++ 
 LINK	= g++
 CXXFLAGS	= -Wcast-align -Waggregate-return -Wall -O4 -funroll-loops
#-DSOCKETIO
 CXXFLAGS_DBUG = -g 
 DEPFLAGS =	-MM
 GXX_LD=/usr/bin/ld
 LIB	= -lm 
 PURIFY =
 QUANTIFY =

ifeq ($(OSTYPE), IRIX)
  LIB	= -lm -lmalloc
endif

ifeq ($(OSTYPE), Solaris2)
  # Use Sparcworks
 CXX	= CC
 LINK	= CC
 DEPFLAGS	= -xM1
 CXXFLAGS	= -ptr$(TDBDIR) -fast 
 CXXFLAGS_DBUG	= -ptr$(TDBDIR) -g -xsb 
#-DSOCKETIO
#  LIB	= -lm -lsocket -lnsl
  PURIFY =	purify -always-use-cache-dir=yes -cache_dir=/tmp/pure_cache
  QUANTIFY =	quantify -always-use-cache-dir=yes -cache_dir=/tmp/quant_cache
endif

CFLAGS = -D$(OSTYPE) -DTIMING 
#-DNDEBUG
#-DLM_STATS

CPPFLAGS = $(CFLAGS) $(CXXFLAGS)  
CPPFLAGS_DBUG =  $(CFLAGS) $(CXXFLAGS_DBUG)


SRCS	=  NW_acoustic.cc NW_scheduler.cc NW_statics.cc NW_decoder.cc \
           NW_subword.cc NW_node.cc \
           NW_acoustic_lexicon.cc NW_tree_lexicon.cc NW_linear_lexicon.cc \
	   NW_hypothesis.cc NW_param.cc NW_misc.cc \
	   NW_lm.cc NW_ngram.cc  NW_tagged_ngram.cc \
	   NW_mixture_ngram.cc \
	   NW_fsn.cc \
	   NW_lattice.cc 
# malloc.cc 

# socketclass.h 
HDRS	=  all.h NW_acoustic.h NW_scheduler.h \
	   NW_hmm.h NW_subword.h NW_node.h \
	   NW_tree_lexicon.h  NW_linear_lexicon.h \
           NW_hypothesis.h NW_stack.h \
           NW_lattice.h  \
	   NW_queue.h NW_type.h NW_pqueue.h \
           NW_hash.h NW_list.h NW_collection.h \
	   NW_lexicon.h NW_acoustic_lexicon.h \
	   NW_lm.h NW_ngram.h NW_fsn.h NW_tagged_ngram.h \
	   NW_mixture_ngram.h \
           NW_misc.h NW_param.h NW_debug.h



OBJS    = $(patsubst %.cc,$(OBJDIR)/%.o,$(SRCS))

# C++ files
$(OBJDIR)/%.o:		%.cc
		$(CXX) $(CPPFLAGS) -c $< -o $@
#		$(CXX) $(CPPFLAGS_DBUG) -c $< -o $@

#all:		noway

noway:		$(OBJS) 
		$(LINK) $(CPPFLAGS) -o $(OBJDIR)/noway $^ $(LIB)
		rm -f noway; ln -s $(OBJDIR)/noway noway

noway.pure:	$(OBJS) 
		$(PURIFY) $(LINK) $(CPPFLAGS_DBUG) -o $(OBJDIR)/noway.pure $^ $(LIB)
		rm -f noway.pure; ln -s $(OBJDIR)/noway.pure noway.pure

noway.quant:	$(OBJS) 
		$(QUANTIFY) $(LINK) $(CPPFLAGS_DBUG) -o $(OBJDIR)/noway.quant $^ $(LIB)
		rm -f noway.quant; ln -s $(OBJDIR)/noway.quant noway.quant


noway.debug:	$(OBJS) 
		$(LINK) $(CPPFLAGS_DBUG) -o $(OBJDIR)/noway.debug $^ $(LIB)
		rm -f noway.debug; ln -s $(OBJDIR)/noway.debug noway.debug

$(OBJS):	Makefile

install:	$(OBJDIR)/noway
		install -m 0775 $< $(BINDIR)

clean:	        
		rm -f $(OBJDIR)/*.o $(OBJDIR)/noway $(OBJDIR)/noway.pure  \
		  $(OBJDIR)/noway.quant

spotless:	clean
		rm -f *.~*~


# Distribution
DIST	=  $(SRCS) $(HDRS) Makefile noway.1 README CHANGES INSTALL OSType BINType 


dist:	        #progman.ps
		mkdir noway-$(VERSION); cp $(DIST) noway-$(VERSION);\
                tar cf - noway-$(VERSION) | gzip > ./noway-$(VERSION).tar.gz;\
		rm -rf noway-$(VERSION)

# Dependency generation
DEPS	= $(patsubst %.cc,$(DEPDIR)/%.d,$(SRCS))

$(DEPDIR)/%.d:		%.cc
		$(SHELL) -ec '$(CXX) $(CFLAGS) $(DEPFLAGS) $< | sed '\''s/$*.o/$(OBJDIR)\/& $(patsubst $(DEPDIR)/%,$(DEPDIR)\/%,$@) /g'\'' > $@'

-include $(DEPS)
