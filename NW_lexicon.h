// -*- Mode: C++;  -*-
// File: lexicon.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Tree-structured lexicon.  Contains phoneset and interface 
//*             to acoustics.
//*
//* CLASSES:   LEXICON
//* 
//* REQUIRES:  HMM (hmm.h);  LIST<T> (list.h); STR_HASH_TABLE<T> (hash.h)
//*
//* HISTORY:
//*  Jan 27 12:58 1995 (sjr): Added support for state path decoding
//*  May 24 15:47 1994 (sjr): Works okay in decoder
//*  Apr 12 14:09 1994 (sjr): Tree structure in place and tested.
//* Created: Wed Apr  6 16:12:40 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef NW_LEXICON_H
#define NW_LEXICON_H

// Tree-structured lexicon, also includes Viterbi routines (to
//  extend a hypothesis)
// All the work is done by this class or it's component nodes.

#include "NW_type.h"
#include "NW_hash.h"
#include "NW_list.h"


class LEXICON {
public:
  LEXICON () {
    vocabulary = new LIST<STRING>(25000);
    unigrams = new LIST<INT>(25000);
    npron = 0;
    nmissing = 0;
    vocab_hash = new STR_HASH_TABLE<INT>(-1);
  }

  virtual ~LEXICON() { 
    delete vocabulary;
    delete unigrams;
    delete vocab_hash;
  }

  // Dictionary and phoneset
  LIST<STRING> *vocabulary;
  LIST<INT> *unigrams;
  STR_HASH_TABLE<INT> *vocab_hash;
  INT npron;
  INT nmissing;

  virtual  VOID eat_acoustic_input(INT) 
    { panic("LEXICON::eat_acoustic_input() not implemented!"); }

};

#endif
