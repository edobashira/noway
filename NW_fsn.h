// -*- Mode: C++;  -*-
// File: NW_fsn.h
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


#include "NW_lm.h"

class FSN_LINK;

class FSN_NODE {
  LIST<FSN_LINK*>* links;
  INT nin;
public:
  FSN_NODE() {
    links = new LIST<FSN_LINK*>();
    nin = 0;
  }

  ~FSN_NODE() {
    delete links;
  }

  INT get_nin() { return nin; }
  VOID inc_nin() { nin++; }
  VOID dec_nin() { nin--; }
  LIST<FSN_LINK*>* get_links() {return(links); }
};


class FSN_LINK {
  INT word;
  INT prob;
  FSN_NODE* start_node;
  FSN_NODE* end_node;
public:
  FSN_LINK(INT w, FSN_NODE* start, FSN_NODE* end, INT p = 0) {
    word = w;
    start_node = start;
    end_node = end;
    prob = p;
  }

  ~FSN_LINK() {
    if(end_node->get_nin() == 1)
      delete end_node;
    else
      end_node->dec_nin();
  }

  INT get_word() { return(word); }
  INT get_prob() { return(prob); }
  FSN_NODE* get_end_node() { return(end_node); }
  FSN_NODE* get_start_node() { return(start_node); }
};



class FSN : public GRAMMAR {
protected:
  FSN_NODE* root;
  LIST<INT>* rev_hist;

  FSN_NODE* current_node(LIST<INT> *history, INT = -1) {
    rev_hist->fast_clear();
    history->init_cursor();
    while(history->next())
      rev_hist->push(history->read());
    return current_node_reverse_history(rev_hist);
  }


  FSN_NODE* current_node(HYP *h, INT = -1) {
    HYP* hyp;

    rev_hist->fast_clear();
    hyp = h;
    while(hyp != NULL) {
      rev_hist->push(hyp->get_word());
      hyp = hyp->get_previous();
    }
    return current_node_reverse_history(rev_hist);
  }

  FSN_NODE* current_node_reverse_history(LIST<INT> *history){
    LIST<FSN_LINK*>* l;
    FSN_NODE* nd = root;
    INT w;
    BOOL fail;

    history->init_cursor();
    while(history->next()) {
      w = history->read();
      fail = TRUE;
      if(nd == NULL)
	break;
      l = nd->get_links();
      l->init_cursor();
      while(l->next()){
	if(l->read()->get_word() == w){
	  nd = l->read()->get_end_node();
	  fail = FALSE;
	  break;
	}
      }
      if(fail == TRUE)
	return(NULL);
    }
    return(nd);
  }

    
  
public:
  FSN() : GRAMMAR() {
    root = new FSN_NODE();
    rev_hist = new LIST<INT>();
  }

  ~FSN() {
    delete root;
  }

  INT build(GENFILE *f, LEXICON *l, INT sc, INT os, BOOL c, 
	    INT nhy, BOOL build_vocab = TRUE);
  INT build_from_ref(QUEUE<INT>* wq);

  // return PROB_NULL if not present
  INT get_prob(LIST<INT> *history, INT nextword, INT hypid = -1) {
    FSN_NODE* nd = current_node(history, hypid);
    return get_prob(nd, nextword);
  }
    
  INT get_prob(HYP* h, INT nextword, INT hypid = -1){
    FSN_NODE* nd = current_node(h, hypid);
    return get_prob(nd, nextword);
  }

  INT get_prob(FSN_NODE* nd, INT nextword) {
    if(nd == NULL)
      return(PROB_NULL);

    LIST<FSN_LINK*>* l;
    l = nd->get_links();
    l->init_cursor();
    while(l->next()) {
      if(l->read()->get_word() == nextword)
	return l->read()->get_prob();
    }
    return PROB_NULL;
  }

  

  LIST<INT>* get_successors(HYP* h);

  // VOID reset()
  // VOID update_cache(LIST<HYP*>* hlist);

  // BOOL hyp_similar(LIST<INT> *h1h, LIST<INT> *h2h);
  // BOOL  hyp_similar(HYP *h1, HYP *h2);
  // BOOL hyp_similar(LIST<INT> *h1h, LIST<INT> *h2h, INT w);
  // BOOL  hyp_similar(HYP *h1, HYP *h2, INT w);

  VOID print(GENFILE *, LEXICON *) {}
  VOID write(GENFILE *, LEXICON *) {}
  VOID info(GENFILE *) {}

  //  DOUBLE text_decode(GENFILE *f, GENFILE *op, LEXICON *l, DOUBLE &noUNKp, 
  //		     BOOL ignore_nl = FALSE);
  //  DOUBLE text_score(QUEUE<INT> *wq, GENFILE *f, LEXICON *l, DOUBLE &noUNK_p, 
  //		   INT &noUNK_l);
};



class ALIGN_LM : public GRAMMAR {
  const INT ALIGN_BUFFER_SIZE;
  QUEUE<FSN*>* utt_models;
  CHAR* albuffer;
public:
  ALIGN_LM () :  ALIGN_BUFFER_SIZE(1024*1024+1) {
    utt_models = new QUEUE<FSN*>();
    albuffer = new CHAR[ALIGN_BUFFER_SIZE];
  }

  ~ALIGN_LM () {
    delete utt_models;
    delete[] albuffer;
  }

  VOID end_utt() { delete utt_models->dequeue(); }

  INT n_utt() { return utt_models->size; }

  INT build(GENFILE *f, LEXICON *l, INT sc = 0, INT os = 0, BOOL c=FALSE, 
	    INT nhy=-1, BOOL build_vocab=TRUE);

  INT get_prob(LIST<INT> *history, INT nextword, INT hyp = -1){ 
    assert(get_utt_models()->size > 0);
    return get_utt_models()->top()->get_prob(history, nextword, hyp); 
  }
  INT get_prob(HYP* h, INT nextword, INT hyp = -1){ 
    assert(get_utt_models()->size > 0);
    return get_utt_models()->top()->get_prob(h, nextword, hyp); 
  }

  VOID info(GENFILE *f)
    { f->fprintf("ALIGN_LM::info() not implemented"); }

  VOID print(GENFILE *, LEXICON *) {}
  VOID print_unigrams(GENFILE *, LEXICON *) {}
  VOID write(GENFILE *, LEXICON *) {}
  DOUBLE text_decode(GENFILE *, GENFILE *, LEXICON *, DOUBLE &noUNKp, 
		     BOOL){
    noUNKp = 0.0;
    return(0.0);
  }

  DOUBLE text_score(QUEUE<INT> *, GENFILE *, LEXICON *, DOUBLE &, INT &)
    {panic("ALIGN_LM::text_score() not implemented!"); return(-1);}

  BOOL  hyp_similar(LIST<INT> *h1, LIST<INT> *h2){ 
    assert(get_utt_models()->size > 0);
    return get_utt_models()->top()->hyp_similar(h1, h2);
  }

  BOOL  hyp_similar(LIST<INT> *h1, LIST<INT> *h2, INT w){ 
    assert(get_utt_models()->size > 0);
    return get_utt_models()->top()->hyp_similar(h1, h2, w);
  }

  BOOL  hyp_similar(HYP *h1, HYP *h2){ 
    assert(get_utt_models()->size > 0);
    return get_utt_models()->top()->hyp_similar(h1, h2);
  }

  BOOL  hyp_similar(HYP *h1, HYP *h2, INT w){ 
    assert(get_utt_models()->size > 0);
    return get_utt_models()->top()->hyp_similar(h1, h2, w);
  }

  QUEUE<FSN*>* get_utt_models() { return utt_models; }

 LIST<INT>* get_successors(HYP* h)
    { return get_utt_models()->top()->get_successors(h); }
};

