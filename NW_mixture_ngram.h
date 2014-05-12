// -*- Mode: C++;  -*-
// File: NW_mixture_ngram.h
// Author: Steve Renals (s.renals@dcs.shef.ac.uk)
// Copyright (C) University of Sheffield, 1997-98

#ifndef NW_MIXTURE_NGRAM_H
#define NW_MIXTURE_NGRAM_H

#include "NW_ngram.h"

class MIXTURE_NGRAM : public NGRAM {
protected:
  LIST<NGRAM*>* ng_lst;
  //  LIST<INT>* ng_probs;

  INT uniform_wt;

  INT nunigrams;

  VOID get_component_probs(HYP* h, INT nextword, INT hypid = -1);

public:
  MIXTURE_NGRAM () : NGRAM() {
    ng_lst = new LIST<NGRAM*>();
    //    ng_probs = new LIST<INT>();
    nunigrams = -1;
    uniform_wt = 0;
  }

  ~MIXTURE_NGRAM () {
    delete ng_lst;
    //    delete ng_probs;
  }

  INT build(GENFILE *f, LEXICON *l, INT sc, INT os, BOOL cache_lm, 
	    INT nhy = -1, BOOL build_vocab = TRUE);

  INT get_nunigrams() { return nunigrams; }
  VOID set_nunigrams(INT n) { nunigrams = n; }
  INT get_ncomponents() { return ng_lst->size; }

  VOID update_cache(LIST<HYP*>* hlist);

  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1){
    panic("only get_prob(HYP*, INT, INT) works for mixture ngrams!\n");
    return(0);
  }


  INT get_prob(HYP* h, INT nextword, INT hypid = -1);


  DOUBLE text_score(QUEUE<INT> *wq, GENFILE *outfp, LEXICON *l, 
		    DOUBLE &noUNK_p, INT &noUNK_l);

  VOID write_header(GENFILE *f, LEXICON *l) {}
  VOID write(GENFILE *f, LEXICON *l) {}
  VOID print(GENFILE *f, LEXICON *l) {}

  VOID print_unigrams(GENFILE *f, LEXICON *l) {}

  VOID info(GENFILE* f=gstdout) {}
};

class BLIND_MIXTURE_NGRAM : public MIXTURE_NGRAM {
public:
  BLIND_MIXTURE_NGRAM () : MIXTURE_NGRAM() { }

  ~BLIND_MIXTURE_NGRAM () { }

  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1){
    panic("only get_prob(HYP*, INT, INT) works for mixture ngrams!\n");
    return(0);
  }

  INT get_prob(HYP* h, INT nextword, INT hypid = -1);

};

class LSA_MIXTURE_NGRAM : public MIXTURE_NGRAM {
public:
  LSA_MIXTURE_NGRAM () : MIXTURE_NGRAM() { }

  ~LSA_MIXTURE_NGRAM () { }

  INT build(GENFILE *f, LEXICON *l, INT sc, INT os, BOOL cache_lm, 
	    INT nhy = -1, BOOL build_vocab = TRUE);


  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1){
    panic("only get_prob(HYP*, INT, INT) works for mixture ngrams!\n");
    return(0);
  }


  INT get_prob(HYP* h, INT nextword, INT hypid = -1);

};


#endif

