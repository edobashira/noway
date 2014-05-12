// -*- Mode: C++;  -*-
// File: grammar.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION: Language models
//*
//* CLASSES:  GRAMMAR, NGRAM, BIGRAM, TRIGRAM, WORDPAIR
//* 
//* REQUIRES: LEXICON (lexicon.h)
//*
//* HISTORY:
//*  Jul 14 12:08 1994 (sjr): BIGRAM and TRIGRAM test okay, with 
//*                           no, full or incremental caching
//*  May 24 15:41 1994 (sjr): Added WORDPAIR class untested
//*  May  4 11:42 1994 (sjr): BIGRAM tested okay.  Read in 5K bigram 
//*                            grammar (~10Mb for 825K bigrams), and
//*                            tested recall and caching.
//* Created: Fri Apr 29 16:13:10 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef NW_LM_H
#define NW_LM_H
// To facilitate storage as USHORTs, all LM probabilities are 
//  stored as -log10(prob).

// This defines the interface for derived grammar
class LEXICON;
class HYP;
class STACK;


class GRAMMAR {
protected:
//  VOID swap2(VOID* ptr){
//    char* c = (char*)ptr;
//    char foo = c[0];
//    c[0] = c[1];
//    c[1] = foo;
//  }

//  VOID swap4(VOID* ptr){
//    char* c = (char*)ptr;
//    char foo = c[0];
//    c[0] = c[3];
//    c[3] = foo;
//    foo = c[1];
//    c[1] = c[2];
//    c[2] = foo;    
//  }

public:
  GRAMMAR() { 
    scale = 1; 
    offset = 0;
    //    nunigrams =0;
    use_unk = FALSE;
    unk_factor = 0;
    tagged = FALSE;
  }

  virtual ~GRAMMAR() {}

  INT scale;
  INT offset;
  //  INT nunigrams;
  BOOL tagged;
  BOOL use_unk;
  USHORT  unk_factor;

  // All grammars must define these
  virtual INT build(GENFILE* f, LEXICON* l, INT sc, INT os, 
		    BOOL cache_lm, INT nhy, BOOL build_vocab = TRUE) = 0;

  virtual INT get_nunigrams() { return 0; }

  // returns log10(prob)
  virtual INT get_prob(LIST<INT> *history, INT nextword, INT hyp = -1) = 0;
  virtual INT get_prob(HYP* h, INT nextword, INT hyp = -1) = 0;

  virtual LIST<INT>* get_successors(HYP*) {
    panic("Not yet implemented!");
    return NULL;
  }

  // redefined by cache classes to clear the cache
  virtual VOID reset(){}
  virtual VOID update_cache(LIST<HYP*>* ){}

  virtual VOID info(GENFILE *) = 0;

  virtual VOID print(GENFILE *, LEXICON *) = 0;
  virtual VOID print_unigrams(GENFILE *, LEXICON *) {}
  virtual VOID write(GENFILE *, LEXICON *) = 0;
  virtual DOUBLE text_decode(GENFILE *f, GENFILE *op, LEXICON *l, DOUBLE &noUNKp, BOOL ignore_nl = FALSE);
  virtual DOUBLE text_score(QUEUE<INT> *, GENFILE *, LEXICON *, DOUBLE &, INT &)
    {panic("text_score() not implemented!"); return(-1);}

  INT end_of_sentence_prob(LIST<INT> *h){
    if(h->top() == Sentence_end_index)
      return(0);
    else
      return(get_prob(h, Sentence_end_index));
  }

  INT end_of_sentence_prob(HYP *h){
    if(h->get_word() == Sentence_end_index)
      return(0);
    else
      return(get_prob(h, Sentence_end_index));
  }

  INT end_of_sentence_prob(INT* h){
    assert(get_order() > 0);
    if(h[0] == Sentence_end_index)
      return(0);
    LIST<INT>* lh = new LIST<INT>();
    for(INT i = get_order()-2; i >= 0; i++)
      lh->push(h[i]);
    return(get_prob(lh, Sentence_end_index));
    delete lh;
  }


  virtual BOOL  hyp_similar(LIST<INT> *h1, LIST<INT> *h2)
    { return hyp_exact_similar(h1, h2); }

  virtual BOOL  hyp_similar(LIST<INT> *h1, LIST<INT> *h2, INT w)
    { return hyp_exact_similar(h1, h2, w); }

  virtual BOOL  hyp_similar(HYP *h1, HYP *h2)
    { return hyp_exact_similar(h1, h2); }
    
  virtual BOOL  hyp_similar(HYP *h1, HYP *h2, INT w)
    { return hyp_exact_similar(h1, h2, w); }

  BOOL hyp_exact_similar(LIST<INT> *h1h, LIST<INT> *h2h);
  BOOL hyp_exact_similar(LIST<INT> *h1h, LIST<INT> *h2h, INT w);

  BOOL hyp_exact_similar(HYP *h1h, HYP *h2h);
  BOOL hyp_exact_similar(HYP *h1h, HYP *h2h, INT w);

  virtual INT oov2unk(INT wdid){
    if(wdid & LM_OOV_MASK) 
      return (Unknown_word_index);
    else 
      return (wdid);
  }

  virtual BOOL is_ngram() { return FALSE; }
  virtual INT get_order() { return 0; }
  virtual INT get_ncomponents() { return 1; }

  virtual INT align(INT, LIST<HYP*>*, LEXICON*, QUEUE<STACK*>*) { return 0; }
};

class NOGRAM : public GRAMMAR {
public:
  NOGRAM() : GRAMMAR() {}

  // All grammars must define these
  INT build(GENFILE*, LEXICON*, INT, INT, BOOL, INT, BOOL) { return(0); }

  // returns log10(prob)
  INT get_prob(LIST<INT>*, INT, INT)
    { return(0); }
  INT get_prob(HYP*, INT, INT)
    { return(0); }



  VOID info(GENFILE *f)
    { f->fprintf("No Language Model\n"); }

  VOID print(GENFILE *, LEXICON *) {}
  VOID print_unigrams(GENFILE *, LEXICON *) {}
  VOID write(GENFILE *, LEXICON *) {}
  DOUBLE text_decode(GENFILE*, GENFILE*, LEXICON*, DOUBLE &noUNKp, 
		     BOOL){
    noUNKp = 0.0;
    return(0.0);
  }

  DOUBLE text_score(QUEUE<INT>*, GENFILE *, LEXICON *, DOUBLE &, INT &)
    {panic("text_score() not implemented!"); return(-1);}

  BOOL hyp_similar(LIST<INT>*, LIST<INT>*)
    { return(TRUE); }
  BOOL  hyp_similar(HYP*, HYP*) 
    { return(TRUE); }

  BOOL hyp_similar(LIST<INT>*, LIST<INT>*, INT)
    { return(TRUE); }
  BOOL  hyp_similar(HYP*, HYP*, INT)
    { return(TRUE); }
    
};
#endif
