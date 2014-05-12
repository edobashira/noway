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

#ifndef NW_ACOUSTIC_LEXICON_H
#define NW_ACOUSTIC_LEXICON_H

#include "NW_lexicon.h"

class ACOUSTIC_LEXICON : public LEXICON {
protected:
  HMM* pause_model;
  INT nlinear_nodes;
  INT silence_idx;
public:

  ACOUSTIC_LEXICON () : LEXICON() {
    phoneset = new STR_HASH_TABLE<HMM*>();
    phonenames = new LIST<STRING>();
    node_count = 0;
    nlinear_nodes = 0;
#ifdef RNN_SILENCE_INDEX
  silence_idx = RNN_SILENCE_INDEX;
#else
  silence_idx = -1;
#endif
  }

  // phoneset
  STR_HASH_TABLE<HMM*> *phoneset;
  LIST<STRING> *phonenames;
  INT node_count;

  LIST<FLOAT> *make_phoneset(PARAM_TABLE *parm, BOOL phi, GENFILE *phonefp, GENFILE *priorfp = NULL); 
  VOID read_phone_models(PARAM_TABLE *parm, GENFILE *phfp, 
			 FLOAT dscale = 1.0, FLOAT pdp = 1.0);
  LIST<FLOAT> *read_priors(GENFILE *prfp);
  LIST<FLOAT> *read_phi(PARAM_TABLE *parm, GENFILE *phfp, 
			FLOAT dscale = 1.0, FLOAT pdp = 1.0);

  // Acoustics
  VOID set_acoustic_input(PARAM_TABLE *p, LIST<FLOAT> *pri);
  VOID eat_acoustic_input(INT t) {
    INT x;
    while(t < pause_model->get_sentence_length()){
      x = pause_model->output_prob(t);
      t++;
    }
  }
  VOID init_lub(INT t) { pause_model->hmm_set_lub(t); }
  INT sentence_length() { return(pause_model->get_sentence_length()); }
  BOOL skip_frames(INT n) { return(pause_model->skip_frames_in_current_sent(n)); }

  // virtual functions
  virtual VOID build(PARAM_TABLE *parm, const INT nw, BOOL phi, GENFILE *lfp, 
		     GENFILE *phfp, GENFILE *prfp = NULL) = 0;
  virtual VOID score(HYP* h, GRAMMAR *lm, LIST<INT> *ps, LIST<INT> *ps2=NULL) = 0;
  virtual INT extend(INT t, LIST<HYP*> *l, GRAMMAR *lm, QUEUE<STACK*> *stks) = 0;
  virtual INT align(INT t, LIST<HYP*> *l, GRAMMAR *lm, QUEUE<STACK*> *stks) = 0;

 };


#endif
