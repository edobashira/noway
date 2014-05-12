// -*- Mode: C++;  -*-
// File: ngram.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* Function: Language models
//*
//* CLASSES:  NGRAM, NGRAM_FINAL_OBJ, NGRAM_INTER_OBJ
//* 
//* REQUIRES: LEXICON (lexicon.h) 
//*
//* HISTORY: Based on noway v1 grammar.h  LM caching is now subclassed
//*  from un
//*
//* Created: Fri May 31 16:13:10 1997 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef NW_TAGGED_NGRAM_H
#define NW_TAGGED_NGRAM_H

#include "NW_ngram.h"
#include "NW_mixture_ngram.h"

class TAG_UNIGRAM {
  // list of arrays of p(w|P)
  //  each array has total_vocab_size elements - all words in pron dictionary, 
  //   including OOvs
  LIST<LIST<INT>*>*  tag_probs;

  // LIST of tag ids (corresponding to LM index)
  LIST<INT>* tag_ids;

  // LIST of tag descriptionms (e.g. LOCATION, UNK, NAME, etc.)
  LIST<CHAR *>* tag_strs;
  
public:
  TAG_UNIGRAM() { 
    tag_probs = new LIST<LIST<INT>*>(); 
    tag_strs = new LIST<char *>(); 
    tag_ids = new LIST<INT>(); 
  }

  ~TAG_UNIGRAM() {
    while(tag_probs->size > 0)
      delete tag_probs->pop();
    delete tag_probs;

    delete tag_ids;

    while(tag_strs->size > 0)
      delete tag_strs->pop();
    delete tag_strs;
  }

  INT build(GENFILE *f, LEXICON *l);

  // returns LM tag inbdex t s.t. ARGMAX_t P(w|t)  
  INT max_prob_tag(INT w) {
    INT indx = w & WORD_MASK;
    INT res = Unknown_word_index;
    INT maxp = -INT_MAX;
    tag_probs->init_cursor();
    while(tag_probs->next()) {
      if (indx < tag_probs->read()->size &&
	  tag_probs->read()->get(indx) != PROB_NULL &&
	  tag_probs->read()->get(indx) > maxp) {
	maxp = tag_probs->read()->get(indx);
	res = tag_ids->get(tag_probs->get_cursor());
      }
    }
    return res;
  }

  VOID prob_init_cursor() {
    tag_probs->init_cursor();
  }

  BOOL prob_next() {
    return tag_probs->next();
  }


  LIST<INT>* prob_read() {
    return(tag_probs->read());
  }

  VOID id_init_cursor() {
    tag_ids->init_cursor();
  }

  BOOL id_next() {
    return tag_ids->next();
  }


  INT id_read() {
    return(tag_ids->read());
  }

  VOID strs_init_cursor() {
    tag_strs->init_cursor();
  }

  BOOL strs_next() {
    return tag_strs->next();
  }


  CHAR* strs_read() {
    return(tag_strs->read());
  }

  // returns p(w | tag) : tag is the LM index for the tag
  INT get_tagged_prob(INT w, INT lmtag) { 
    INT indx = get_tag_index(lmtag);
    if(indx == -1)
      panic("TAG_UNIGRAM::get_tagged_prob() - Unknown tag index (%d)\n", lmtag);
    return(tag_probs->get(indx)->get(w&WORD_MASK)); 
  }  

  STRING get_tag_str(INT t) 
    { return(tag_strs->get(t)); }

  // returns TAG index of tagstr
  INT get_tag_index(STRING tagstr) {
    tag_strs->init_cursor();
    while(tag_strs->next()) {
      if(!strcmp(tag_strs->read(), tagstr))
	tag_strs->get_cursor();
    }
    return(-1);
  }

  // returns  TAG index of lmtag
  INT get_tag_index(INT lmtag) {
    tag_ids->init_cursor();
    while(tag_ids->next()) {
      if(tag_ids->read() == lmtag)
	return tag_ids->get_cursor();
    }
    return(-1);
  }

  VOID info(GENFILE *f);

};

class TAGGED_NGRAM : public GRAMMAR {
protected:
  NGRAM *ng;
  TAG_UNIGRAM *tag;

public:
  TAGGED_NGRAM () : GRAMMAR() {
    ng = new NGRAM();
    tag = new TAG_UNIGRAM();
    tagged = TRUE;
  }

  ~TAGGED_NGRAM () {
    delete ng;
    delete tag;
  }

  INT build(GENFILE *f, LEXICON *l, INT sc, INT os, BOOL cache_lm, INT
	    nhy = -1, BOOL build_vocab = TRUE){
    if(cache_lm == TRUE) {
      fprintf(stderr, "LM caching not implemented for TAGGED_NGRAM\n");
      cache_lm = FALSE;
    }
    return(ng->build(f, l, sc, os, cache_lm, nhy, build_vocab));
  }

  INT build_tags(GENFILE *f, LEXICON *l) {
    return(tag->build(f, l));
  }

  INT get_order() { return ng->get_order(); }

  VOID write_header(GENFILE *f, LEXICON *l) {}
  VOID write(GENFILE *f, LEXICON *l) {}
  VOID print(GENFILE *f, LEXICON *l) {}

  VOID print_unigrams(GENFILE *f, LEXICON *l) {}

  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1);
  INT get_prob(HYP* h, INT nextword, INT hypid = -1);

  DOUBLE text_score(QUEUE<INT> *wq, GENFILE *outfp, LEXICON *l, 
  		 DOUBLE &noUNK_p, INT &noUNK_l);
    
  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2) {
    return ng->hyp_similar(h1, h2);
  }
  BOOL hyp_similar(HYP *h1, HYP *h2){
    return ng->hyp_similar(h1, h2);
  }

  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2, INT w){
    h2->push(oov2unk(w));
    BOOL res = hyp_similar(h1,h2);
    h2->pop();
    return(res);    
  }

  BOOL hyp_similar(HYP *h1, HYP *h2, INT w){
    return ng->hyp_similar(h1, h2, oov2unk(w));
  }

  VOID reset(){ ng->reset(); }
  VOID update_cache(LIST<HYP*>* hl) { ng->update_cache(hl); }

  VOID info(GENFILE* f=gstdout) { ng->info(f);  tag->info(f); }

  INT oov2unk(INT w) {
    if(w & LM_OOV_MASK) {
      return(tag->max_prob_tag(w));
    } else {
      return (w);
    }
  }

};


class TAGGED_MIXTURE_NGRAM : public TAGGED_NGRAM {
public:
  TAGGED_MIXTURE_NGRAM () {
    scale = 1; 
    offset = 0;
    use_unk = FALSE;
    unk_factor = 0;
    tagged = TRUE;
    ng = new MIXTURE_NGRAM();
    tag = new TAG_UNIGRAM();
  }

  ~TAGGED_MIXTURE_NGRAM () {
    delete ng;
    delete tag;
  }

};


class TAGGED_BLIND_MIXTURE_NGRAM : public TAGGED_MIXTURE_NGRAM {
public:
  TAGGED_BLIND_MIXTURE_NGRAM () {
    scale = 1; 
    offset = 0;
    use_unk = FALSE;
    unk_factor = 0;
    tagged = TRUE;
    ng = new BLIND_MIXTURE_NGRAM();
    tag = new TAG_UNIGRAM();
  }

  ~TAGGED_BLIND_MIXTURE_NGRAM () {
    delete ng;
    delete tag;
  }

};


#endif
