
// -*- Mode: C++;  -*-
// File: all.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Header files for noway
//*
//* CLASSES:
//* 
//* RELATED PACKAGES:
//*
//* HISTORY:
//* Created: Tues Apr 5 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// switch off assertions
// #define NDEBUG
#include <assert.h>
#include <ctype.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
//#include <sys/mman.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TIMING
#include <time.h>
#include <sys/times.h>
#endif
#include <unistd.h>
#if 0
#if defined(Solaris2) || defined(Linux)
#include <sys/time.h>
#include <sys/resource.h>
#endif/* Solaris2 Linux */
#if defined(Solaris2)
extern int getrusage(int, struct rusage *);
#endif /* Solaris2 */
#endif

#ifdef SOCKETIO
#include "socketclass.h"
#endif
#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_param.h"

#include "NW_debug.h"
#include "NW_hypothesis.h"
#include "NW_lm.h"
#include "NW_ngram.h"
#include "NW_mixture_ngram.h"
#include "NW_tagged_ngram.h"
#include "NW_fsn.h"
#include "NW_lattice.h"
#include "NW_stack.h"

// Decoder classes
#include "NW_acoustic.h"
#include "NW_hmm.h"
#include "NW_subword.h"
#include "NW_node.h"
#include "NW_tree_lexicon.h"
#include "NW_linear_lexicon.h"



