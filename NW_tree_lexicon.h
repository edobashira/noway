// -*- Mode: C++;  -*-
// File: tree_lexicon.h
// Author: Steve Renals (s.renals@dcs.shef.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
// Copyright (C) University of Sheffield, 1994-97
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


// Tree-structured lexicon, also includes Viterbi routines (to
//  extend a hypothesis)
// All the work is done by this class or it's component nodes.
class GRAMMAR;

#include "NW_acoustic_lexicon.h"

class TREE_LEXICON : public ACOUSTIC_LEXICON {
  // The nodes are the tree structured lexicon need only be accessed by
  //  the lexicon
protected:
  NODE* root;
  INT nlevels;

  //  QUEUE<NODE*> *active;
  //  QUEUE<NODE*> *next_active;
  //  MINPQUEUE<NODE*> *exceptions;
  QUEUE<QUEUE<NODE*>*>* active;
  QUEUE<QUEUE<NODE*>*>* next_active;

public:
  TREE_LEXICON () : ACOUSTIC_LEXICON() {
    pwcount = 0;
    active = new QUEUE<QUEUE<NODE*>*>();
    next_active = new QUEUE<QUEUE<NODE*>*>();
    //    exceptions = new MINPQUEUE<NODE*>();
    nlevels = 0; // root
    depth_count  = new QUEUE<INT>();
  }


  // Construct
  VOID build(PARAM_TABLE *parm, const INT nw, BOOL phi, GENFILE *lfp, GENFILE *phfp, GENFILE *prfp = NULL); 
  VOID allocate_maxlmprob(INT nh);
  VOID commission();

  // Decode
  INT beam() { return(root->beam); }
  //VOID new_sentence() { root->new_sentence(); }
  VOID init_extend(INT t, LIST<HYP*> *l, GRAMMAR *lm);
  VOID score(HYP* h, GRAMMAR *lm, LIST<INT> *ps, LIST<INT> *ps2=NULL);
  INT extend(INT t, LIST<HYP*> *l, GRAMMAR *lm, QUEUE<STACK*> *stks);
  INT align(INT, LIST<HYP*>*, GRAMMAR*, QUEUE<STACK*>*)
    { panic("TREE_LEXICON::align() not implemented\n"); return 0; }

  // Debugging
  VOID info();
  VOID dump();

  QUEUE<INT> *depth_count;
  INT pwcount;
};




