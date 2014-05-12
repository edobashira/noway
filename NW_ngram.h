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

#ifndef NW_NGRAM_H
#define NW_NGRAM_H

#include "NW_lm.h"

// To facilitate storage as USHORTs, all LM probabilities are 
//  stored as -log10(prob).

class NGRAM_FINAL_OBJ;
class NGRAM_INTER_OBJ;

// This is an abstract base class used for arrays of LM nodes which are either 
// terminal (FINAL) or non-terminal (INTER).  Specialised at the bottom 
// to NGRAM_FINAL_ARRAY and  NGRAM_INTER_ARRAY
class NGRAM_ARRAY {
protected:
  INT size;
public:
  NGRAM_ARRAY() {}
  virtual ~NGRAM_ARRAY() {}

  virtual NGRAM_FINAL_OBJ* get_data() = 0;
  virtual NGRAM_FINAL_OBJ* get(INT) = 0;
  virtual VOID load(NGRAM_ARRAY*, INT) = 0;

  VOID resize(INT n){
    if(n > size) {
      panic("NGRAM_ARRAY::resize(%d) --- size = %d, cannot increase size", 
	    n, size);
    }
    size = n;
  }

  INT get_size() { return size; }
};

class TRIGRAM_OBJ {
private:
  USHORT word;
  USHORT prob;
public:
  TRIGRAM_OBJ() {
    word = 0xffff;
    prob = 0;
  }

  VOID init(USHORT w, USHORT p)
    { set_word(w); set_prob(p); }

  INT get_prob() { return(-(INT)prob);  }
  USHORT get_prob_ushort() { return(prob); }
  INT get_word() { return((INT)word); }

  VOID set_word(USHORT w) { word = w; }
  VOID set_prob(USHORT p) { prob = p; }


  USHORT read(GENFILE *infp, BOOL use_unk, USHORT unk_factor);

  VOID write(GENFILE *outfp);
  VOID print(GENFILE *f, LEXICON *l);

  VOID info(GENFILE *f=gstdout) {
    f->fprintf("TRIGRAM_OBJ\n");
    f->fprintf("  Word=%d\n", word);
    f->fprintf("  Prob=%d\n", prob);
  }
};

class BIGRAM_OBJ {
private:
  USHORT word;
  USHORT prob;
  SHORT backoff;
  USHORT n_trigrams;
  TRIGRAM_OBJ* trigrams;
public:
  BIGRAM_OBJ() {
    n_trigrams = 0;
    trigrams = NULL;
  }

  VOID init(USHORT w, USHORT p, SHORT b)
    { set_word(w); set_prob(p); set_backoff(b); }


  INT get_prob() { return(-(INT)prob);  }
  USHORT get_prob_ushort() { return(prob); }
  INT get_word() { return((INT)word); }
  INT get_backoff() { return((INT)backoff); }
  SHORT get_backoff_short() { return(backoff); }

  TRIGRAM_OBJ* get_trigrams() { return(trigrams); }
  INT get_n_trigrams() { return(n_trigrams); }

  VOID set_word(USHORT w) { word = w; }
  VOID set_prob(USHORT p) { prob = p; }
  VOID set_backoff(SHORT b) { backoff = b; }

  TRIGRAM_OBJ* get_trigram_by_index(INT i)
    { assert(i < n_trigrams); return(trigrams+i);}

  TRIGRAM_OBJ* get_trigram(INT w);
  TRIGRAM_OBJ* get_trigram_seq(INT w);
  VOID load_trigrams(INT n, TRIGRAM_OBJ*); 

  VOID clear();

  USHORT read(GENFILE *infp, BOOL use_unk, USHORT unk_factor);

  VOID write(GENFILE *outfp, BOOL nobackoff = FALSE);
  VOID print(GENFILE *f, LEXICON *l, BOOL nobackoff = FALSE);

  VOID info(GENFILE *f=gstdout) {
    f->fprintf("BIGRAM_OBJ\n");
    f->fprintf("  Word=%d\n", word);
    f->fprintf("  Prob=%d\n", prob);
    f->fprintf("  Backoff=%d\n", prob);
    f->fprintf("  %d Trigrams (at %p)\n", n_trigrams, trigrams);
  }  
};

class UNIGRAM_OBJ {
private:
  USHORT word;
  USHORT prob;
  SHORT backoff;
  USHORT n_bigrams;
  BIGRAM_OBJ* bigrams;
public:
  UNIGRAM_OBJ() {
    n_bigrams = 0;
    bigrams = NULL;
  }

  VOID init(USHORT w, USHORT p, SHORT b)
    { set_word(w); set_prob(p); set_backoff(b); }


  INT get_prob() { return(-(INT)prob);  }
  USHORT get_prob_ushort() { return(prob); }
  INT get_word() { return((INT)word); }
  INT get_backoff() { return((INT)backoff); }
  SHORT get_backoff_short() { return(backoff); }

  BIGRAM_OBJ* get_bigrams() { return(bigrams); }
  INT get_n_bigrams() { return(n_bigrams); }

  VOID set_word(USHORT w) { word = w; }
  VOID set_prob(USHORT p) { prob = p; }
  VOID set_backoff(SHORT b) { backoff = b; }

  BIGRAM_OBJ* get_bigram_by_index(INT i)
    { assert(i < n_bigrams); return(bigrams+i);}

  BIGRAM_OBJ* get_bigram(INT w);
  BIGRAM_OBJ* get_bigram_seq(INT w);
  VOID load_bigrams(INT n, BIGRAM_OBJ*);

  VOID clear();

  USHORT read(GENFILE *infp,  BOOL use_unk, USHORT unk_factor);

  VOID write(GENFILE *outfp,  BOOL nobgbackoff = FALSE);
  VOID print(GENFILE *f, LEXICON *l, BOOL recurs = TRUE, BOOL nobgbackoff = FALSE);

  VOID info(GENFILE *f=gstdout) {
    f->fprintf("UNIGRAM_OBJ\n");
    f->fprintf("  Word=%d\n", word);
    f->fprintf("  Prob=%d\n", prob);
    f->fprintf("  Backoff=%d\n", prob);
    f->fprintf("  %d Bigrams (at %p)\n", n_bigrams, bigrams);
  }  
};

class NGRAM : public GRAMMAR {
 protected:

  // Array size nunigrams.  Each bigram[i] contains a list of prob(j) = bigram(i,j) 
  //   and a list bo(j) = backoff(i, j)  and a pointer 
  //  INT nunigrams;
  NGRAM_ARRAY* unigrams;
  INT order;
  INT *context;
  INT *ngram_count;

  INT get_prob_backoff(INT clen);

  CHAR *arpa_hdr_chars;
  //  VOID read_arpa_unigrams(GENFILE *fp, LEXICON *lex, INT sc, INT os, INT ord, 
  //			  INT max_nhyps, BOOL build_vocab = TRUE);

  // Caching

  // cache is max_hyps x nunigrams
  //  cache[i] is array of size nunigrams
  INT** cache;
  INT** tmp_cache;

  // cache_context is max_hyps x order
  INT**  cache_context;
  INT**  tmp_cache_context;

  // backoff cache is max_hyps
  INT*  backoff_cache;
  INT*  tmp_backoff_cache;

  INT cache_size;
  INT tmp_cache_size;

  INT   max_hyps;

  BOOL caching;
  BOOL cached;
public:
  NGRAM() : GRAMMAR() { 
    arpa_hdr_chars = strdup("\\data\\");
    unigrams = NULL;
    context = NULL;
    ngram_count = NULL;
    order = 0;
    //    nunigrams = 0;
    cached = FALSE;
    max_hyps = 0;
    cache_size = 0;
    cache = tmp_cache = NULL;
    cache_context = tmp_cache_context = NULL;
    backoff_cache = tmp_backoff_cache = NULL;
  }

  INT build(GENFILE*f, LEXICON *l, INT sc, INT os, BOOL cache_lm, 
	    INT nhy = -1, BOOL build_vocab = TRUE);

  VOID write_header(GENFILE *f, LEXICON *l);
  VOID write(GENFILE *f, LEXICON *l);
  VOID print(GENFILE *f, LEXICON *l);

  VOID print_unigrams(GENFILE *f, LEXICON *l);

  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1);
  INT get_prob(HYP* h, INT nextword, INT hypid = -1);

  INT get_nunigrams() { return unigrams->get_size(); }

  DOUBLE text_score(QUEUE<INT> *wq, GENFILE *outfp, LEXICON *l, 
  		 DOUBLE &noUNK_p, INT &noUNK_l);
    
  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2);
  BOOL hyp_similar(HYP *h1, HYP *h2);

  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2, INT w);
  BOOL hyp_similar(HYP *h1, HYP *h2, INT w);

  VOID reset() {
    if(cached)
      clear_cache();
  }

  VOID info(GENFILE *outfp=gstdout);

  BOOL is_ngram() { return TRUE; }

  // not defined by LM
  INT read_old_bin_lm(GENFILE *f,  LEXICON *l);
  INT read_arpa_lm(GENFILE *f, LEXICON *l, INT ch, BOOL cache_lm, 
  		   BOOL build_vocab = TRUE);

  VOID update_cache(LIST<HYP*>* hlist);


  VOID create_cache();
  VOID clear_cache();
  INT in_tmp_cache(HYP* h);
  NGRAM_INTER_OBJ* get_context_obj(HYP* h);

  INT get_order() { return order; }

};

class TRIGRAM : public NGRAM {
protected:
  UNIGRAM_OBJ* tg_unigrams;
  INT nunigrams;

public:
  TRIGRAM() : NGRAM() {
    tg_unigrams = NULL;
    order = 3;
    unigrams = 0;
  }

  INT build(GENFILE *f, LEXICON *l, INT sc, INT os, BOOL cache_lm, 
	    INT nhy = -1, BOOL build_vocab = TRUE);
  VOID update_cache(LIST<HYP*>* hlist);

  VOID write(GENFILE*, LEXICON*);
  VOID print(GENFILE*, LEXICON*)
    { panic("use -ngram for printing LMs"); }
  VOID print_unigrams(GENFILE*, LEXICON*)
    { panic("use -ngram for printing LMs"); }

  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1);
  INT get_prob(HYP* h, INT nextword, INT hypid = -1);
  INT get_prob(INT w1, INT w2, INT w3);
  INT get_prob(INT w1, INT w2);

  INT get_nunigrams() { return nunigrams; }
  VOID set_nunigrams(INT n) { nunigrams = n; }

  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2);
  BOOL hyp_similar(HYP *h1, HYP *h2){
    return(h1->get_word() == h2->get_word() &&
	   h1->get_prev_word() == h2->get_prev_word());
  }

  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2, INT w);
  BOOL hyp_similar(HYP *h1, HYP *h2, INT w){
    return(w == h1->get_word() && h1->get_prev_word() == h2->get_word());
  }

  DOUBLE text_score(QUEUE<INT>*, GENFILE*, LEXICON*, 
		    DOUBLE&, INT&) 
    { panic("Use -ngram for text decode!\n"); return 0.0; }

  VOID info(GENFILE* = gstdout) {}

  BIGRAM_OBJ* get_bigram_obj(HYP* h)
    { return tg_unigrams[h->get_prev_word()].get_bigram(h->get_word()); }

};


class BIGRAM : public NGRAM {
protected:
  UNIGRAM_OBJ* bg_unigrams;
  INT nunigrams;

public:
  BIGRAM() : NGRAM() {
    bg_unigrams = NULL;
    order = 2;
    nunigrams = 0;
  }

  INT build(GENFILE *f, LEXICON *l, INT sc, INT os, BOOL cache_lm, 
	    INT nhy = -1, BOOL build_vocab = TRUE);
  VOID update_cache(LIST<HYP*>* hlist);

  VOID write(GENFILE*, LEXICON*);
  VOID print(GENFILE*, LEXICON*)
    { panic("use -ngram for printing LMs"); }
  VOID print_unigrams(GENFILE*, LEXICON*)
    { panic("use -ngram for printing LMs"); }

  INT get_prob(LIST<INT>* history, INT nextword, INT hypid = -1);
  INT get_prob(HYP* h, INT nextword, INT hypid = -1);
  INT get_prob(INT w1, INT w2);

  INT get_nunigrams() { return nunigrams; }
  VOID set_nunigrams(INT n) { nunigrams = n; }

  BOOL hyp_similar(LIST<INT> *h1, LIST<INT> *h2)
    { return(h1->top() == h2->top()); }
  BOOL hyp_similar(HYP *h1, HYP *h2)
    { return(h1->get_word() == h2->get_word()); }

  BOOL hyp_similar(LIST<INT>*h1, LIST<INT>*, INT w)
    { return(w == h1->top()); }
  BOOL hyp_similar(HYP *h1, HYP*, INT w)
    { return(w == h1->get_word() ); }

  DOUBLE text_score(QUEUE<INT>*, GENFILE*, LEXICON *, 
		    DOUBLE&, INT&) 
    { panic("Use -ngram for text decode!\n"); return 0.0; }

  VOID info(GENFILE* = gstdout) {}

};


enum LM_FORMAT {NEWBIN, OLDBIN, ARPA};
enum NGRAM_TYPE {FINAL, INTER};



class NGRAM_FINAL_OBJ {
protected:
  USHORT word;
  USHORT prob;

public:
  NGRAM_FINAL_OBJ() {word = LM_NULL; prob = LM_NULL; }

  virtual NGRAM_TYPE ngram_mode(){ return(FINAL); }

  VOID init(USHORT w, USHORT p, SHORT b);

  INT get_prob() { return(-(INT)prob);  }
  USHORT get_prob_ushort() { return(prob); }
  INT get_word() { return((INT)word); }

  virtual INT get_backoff() { return(0); }
  virtual SHORT get_backoff_short() { return(0); }

  virtual NGRAM_FINAL_OBJ* get_nextgrams() { return(NULL); }
  virtual INT get_n_nextgrams() { return(0); }

  VOID set_word(USHORT w) { word = w; }
  VOID set_prob(USHORT p) { prob = p; }
  virtual VOID set_backoff(SHORT) {}

  virtual NGRAM_FINAL_OBJ* get_nextgram_by_index(INT){ return(NULL); }
  virtual VOID clear() { }

  virtual NGRAM_FINAL_OBJ* get_nextgram(INT){ return(NULL); }
  virtual NGRAM_FINAL_OBJ* get_nextgram_seq(INT){ return(NULL); }
  virtual VOID load_nextgrams(INT , NGRAM_ARRAY*) { }

  virtual USHORT read(GENFILE *infp,  BOOL use_unk, USHORT unk_factor,
		      INT depth_left);

  virtual VOID write(GENFILE *outfp);
  virtual VOID print(GENFILE *f, LEXICON *l, INT level = 0, BOOL recurs=TRUE);

  virtual VOID info(GENFILE *f=gstdout) {
    f->fprintf("NGRAM_FINAL_OBJ\n");
    f->fprintf("  Word=%d\n", word);
    f->fprintf("  Prob=%d\n", prob);
  }
};

class NGRAM_INTER_OBJ : public NGRAM_FINAL_OBJ {
protected:
  SHORT backoff;
  //  USHORT n_nextgrams;
  NGRAM_ARRAY* nextgrams;

public:
  
  NGRAM_INTER_OBJ() : NGRAM_FINAL_OBJ() { 
    backoff = 0; 
    // n_nextgrams = 0; 
    nextgrams = NULL; 
  }

  NGRAM_TYPE ngram_mode(){ return(INTER); }

  INT get_backoff() { return((INT)backoff); }
  SHORT get_backoff_short() { return(backoff); }

  NGRAM_FINAL_OBJ* get_nextgrams() { return(nextgrams->get_data()); }
  INT get_n_nextgrams() { 
    if(nextgrams == NULL) return 0;
    return nextgrams->get_size();
  }

  VOID set_backoff(SHORT b) {backoff = b;}

  NGRAM_FINAL_OBJ* get_nextgram_by_index(INT i)
    { assert(i < get_n_nextgrams()); return(nextgrams->get(i));}

  VOID clear();

  NGRAM_FINAL_OBJ* get_nextgram(INT w);

  // This starts looking from the last bigram access and returns NULL if failure
  // Calling with -1 resets index
  // Assumes ordered list
  NGRAM_FINAL_OBJ* get_nextgram_seq(INT w);

  VOID load_nextgrams(INT n, NGRAM_ARRAY* ng);

  USHORT read(GENFILE *infp,  BOOL use_unk, USHORT unk_factor, 
	      INT depth_left);

  VOID write(GENFILE *outfp);

  VOID print(GENFILE *f, LEXICON *l, INT level = 0, BOOL recurs=FALSE);

  VOID info(GENFILE *f=gstdout) {
    f->fprintf("NGRAM_INTER_OBJ\n");
    f->fprintf("  Word=%d\n", word);
    f->fprintf("  Prob=%d\n", prob);
    f->fprintf("  Backoff=%d\n", backoff);
    f->fprintf("  n_nextgrams = %d\n", get_n_nextgrams());
  }
    

};




class NGRAM_FINAL_ARRAY : public NGRAM_ARRAY {
  NGRAM_FINAL_OBJ* fdata;
public:
  NGRAM_FINAL_ARRAY() { fdata = NULL; }
  NGRAM_FINAL_ARRAY(INT n) { fdata = new NGRAM_FINAL_OBJ[n]; size = n; }
  ~NGRAM_FINAL_ARRAY() { if(fdata != NULL) delete[] fdata;}
  
  NGRAM_FINAL_OBJ* get_data() { return(fdata); }
  NGRAM_FINAL_OBJ* get(INT i) { assert(i < size); return(&fdata[i]); }
  VOID load(NGRAM_ARRAY* ng, INT n){ 
    assert(n <= size);
    memcpy(fdata, ng->get_data(), n*sizeof(NGRAM_FINAL_OBJ)); 
  }
  
};

class NGRAM_INTER_ARRAY : public NGRAM_ARRAY {
  NGRAM_INTER_OBJ* idata;
public:
  NGRAM_INTER_ARRAY() { idata = NULL;}
  NGRAM_INTER_ARRAY(INT n) { idata = new NGRAM_INTER_OBJ[n]; size = n; }
  ~NGRAM_INTER_ARRAY() {if(idata!=NULL) delete[] idata;}

  NGRAM_FINAL_OBJ* get_data() { return(idata); }
  NGRAM_FINAL_OBJ* get(INT i) { assert(i < size); return(&idata[i]); }

  VOID load(NGRAM_ARRAY* ng, INT n){ 
    assert(n <= size);
    memcpy(idata, ng->get_data(), n*sizeof(NGRAM_INTER_OBJ));
  }
};

#endif
