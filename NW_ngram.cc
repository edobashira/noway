// -*- Mode: C++;  -*-
// File: ngram.cc
// Author: Steve Renals (sjr@dcs.shef.ac.uk)
// Copyright (C) Department of Computer Science, Sheffield University, 1997
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* REQUIRES:
//*
//* HISTORY:
//* Created: 
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_debug.h"
#include "NW_param.h"

#include "NW_hypothesis.h"
#include "NW_ngram.h"
#include "NW_lexicon.h"

USHORT
TRIGRAM_OBJ::read(GENFILE *infp, BOOL use_unk, USHORT unk_factor){
  if(infp->fread(&word, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
    //	panic("Ngram %d\n", ncreated);
  }
  if(infp->fread(&prob, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
    //	panic("Ngram %d\n", ncreated);
  }

  if(use_unk && word==Unknown_word_index)
    prob = MIN(65534, unk_factor + prob);
      
  return(prob);
}

VOID
TRIGRAM_OBJ::write(GENFILE *outfp){
  outfp->fwrite(&word, sizeof(USHORT), 1);
  outfp->fwrite(&prob, sizeof(USHORT), 1);
}

VOID
TRIGRAM_OBJ::print(GENFILE *f, LEXICON *l){
  f->fprintf("      ");
  f->fprintf("    %s(%d) %f(%d)n", l->vocabulary->get(word), word, 
	     logdecodelog(get_prob()), get_prob());
}

TRIGRAM_OBJ*
BIGRAM_OBJ::get_trigram(INT w){
  UINT i, lower, upper;
  INT diff;
  
  if(n_trigrams == 0) return(NULL);

  lower = 0;
  upper = n_trigrams;
  i = upper/2;
  while(upper > lower) {
    if((diff = w - trigrams[i].get_word()) == 0)
      return(trigrams+i);
    else if(diff > 0)
      lower = i + 1;
    else // diff < 0
      upper = i;
    i = (lower + upper)/2;
  }
  return(NULL);
}

TRIGRAM_OBJ*
BIGRAM_OBJ::get_trigram_seq(INT w){
  static INT index = 0;
  TRIGRAM_OBJ* res = NULL;
  INT diff;
  
  if(w == -1) { index = 0; return(NULL); }
  if(n_trigrams == 0) return(NULL);

  while(index < n_trigrams){
    if((diff = w - trigrams[index].get_word()) == 0) {
      res = trigrams+index;
      index++; // static var
      break;
    }
    else if(diff > 0) {
      index++;
    }
    else { // diff < 0
      break;
    }
  }
  return(res);

}

VOID
BIGRAM_OBJ::load_trigrams(INT n, TRIGRAM_OBJ* tgo){
  if(n == 0) return;

  assert(trigrams == NULL && n_trigrams == 0);
  n_trigrams = n;
  trigrams = new TRIGRAM_OBJ[n];
  memcpy(trigrams, tgo, n*sizeof(TRIGRAM_OBJ));
}

VOID
BIGRAM_OBJ::clear(){
  if(trigrams != NULL){
    delete trigrams;
    trigrams = NULL;
  }
  n_trigrams = 0;
}

USHORT
BIGRAM_OBJ::read(GENFILE *infp, BOOL use_unk, USHORT unk_factor){
  assert(trigrams == NULL && n_trigrams == 0);
  if(infp->fread(&word, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
  }
  if(infp->fread(&prob, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
  }
  if(infp->fread(&backoff, sizeof(SHORT), 1) < 1){
    perror("Reading trigram");
  }
  if(infp->fread(&n_trigrams, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
  }

  if(use_unk && word==Unknown_word_index)
    prob = MIN(65534, unk_factor + prob);

  if(n_trigrams > 0) {
    trigrams = new TRIGRAM_OBJ[n_trigrams];
    USHORT p;
    for(INT i = 0; i < n_trigrams; i++) 
      p = trigrams[i].read(infp, use_unk, unk_factor);
  }
  return(prob);
}


VOID
BIGRAM_OBJ::write(GENFILE *outfp, BOOL nobackoff){
  outfp->fwrite(&word, sizeof(USHORT), 1);
  outfp->fwrite(&prob, sizeof(USHORT), 1);
  if(!nobackoff){
    outfp->fwrite(&backoff, sizeof(SHORT), 1);
    outfp->fwrite(&n_trigrams, sizeof(USHORT), 1);
  }

  if(!nobackoff)
    for(INT i = 0; i < n_trigrams; i++)
      trigrams[i].write(outfp);
}

VOID
BIGRAM_OBJ::print(GENFILE *f, LEXICON *l, BOOL nobackoff){
  f->fprintf("   ");
  if(nobackoff){
    f->fprintf("    %s(%d) %f(%d)\n", 
	     l->vocabulary->get(word), word, 
	     logdecodelog(get_prob()), get_prob());
  } else {
    f->fprintf("    %s(%d) %f(%d)  %f(%d) (%d trigrams)\n", 
	     l->vocabulary->get(word), word, 
	     logdecodelog(get_prob()), get_prob(),
	     logdecodelog(get_backoff()), get_backoff(),
	     n_trigrams);

    for(INT i = 0; i < n_trigrams; i++)
      trigrams[i].print(f,l);  
  }
}


BIGRAM_OBJ*
UNIGRAM_OBJ::get_bigram(INT w){
  UINT i, lower, upper;
  INT diff;
  
  if(n_bigrams == 0) return(NULL);

  lower = 0;
  upper = n_bigrams;
  i = upper/2;
  while(upper > lower) {
    if((diff = w - bigrams[i].get_word()) == 0)
      return(bigrams+i);
    else if(diff > 0)
      lower = i + 1;
    else // diff < 0
      upper = i;
    i = (lower + upper)/2;
  }
  return(NULL);
}

BIGRAM_OBJ*
UNIGRAM_OBJ::get_bigram_seq(INT w){
  static INT index = 0;
  BIGRAM_OBJ* res = NULL;
  INT diff;
  
  if(w == -1) { index = 0; return(NULL); }
  if(n_bigrams == 0) return(NULL);

  while(index < n_bigrams){
    if((diff = w - bigrams[index].get_word()) == 0) {
      res = bigrams+index;
      index++; // static var
      break;
    }
    else if(diff > 0) {
      index++;
    }
    else { // diff < 0
      break;
    }
  }
  return(res);
}

VOID
UNIGRAM_OBJ::load_bigrams(INT n, BIGRAM_OBJ* bgo){
  if(n == 0) return;

  assert(bigrams == NULL && n_bigrams == 0);
  n_bigrams = n;
  bigrams = new BIGRAM_OBJ[n];
  memcpy(bigrams, bgo, n*sizeof(BIGRAM_OBJ));
}

VOID
UNIGRAM_OBJ::clear(){
  if(bigrams != NULL){
    delete bigrams;
    bigrams = NULL;
  }
  n_bigrams = 0;
}

USHORT
UNIGRAM_OBJ::read(GENFILE *infp, BOOL use_unk, USHORT unk_factor){
  assert(bigrams == NULL && n_bigrams == 0);
  if(infp->fread(&word, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
  }
  if(infp->fread(&prob, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
  }
  if(infp->fread(&backoff, sizeof(SHORT), 1) < 1){
    perror("Reading trigram");
  }
  if(infp->fread(&n_bigrams, sizeof(USHORT), 1) < 1){
    perror("Reading trigram");
  }

  if(use_unk && word==Unknown_word_index)
    prob = MIN(65534, unk_factor + prob);

  if(n_bigrams > 0) {
    bigrams = new BIGRAM_OBJ[n_bigrams];
    USHORT p;
    for(INT i = 0; i < n_bigrams; i++) 
      p = bigrams[i].read(infp, use_unk, unk_factor);
  }
  return(prob);
}

VOID
UNIGRAM_OBJ::write(GENFILE *outfp, BOOL nobgbackoff){
  outfp->fwrite(&word, sizeof(USHORT), 1);
  outfp->fwrite(&prob, sizeof(USHORT), 1);
  outfp->fwrite(&backoff, sizeof(SHORT), 1);
  outfp->fwrite(&n_bigrams, sizeof(USHORT), 1);

  for(INT i = 0; i < n_bigrams; i++)
    bigrams[i].write(outfp,  nobgbackoff);
}

VOID
UNIGRAM_OBJ::print(GENFILE *f, LEXICON *l, BOOL recurs, BOOL noback){
  f->fprintf("    %s(%d) %f(%d)  %f(%d) (%d bigrams)\n", 
	     l->vocabulary->get(word), word, 
	     logdecodelog(get_prob()), get_prob(),
	     logdecodelog(get_backoff()), get_backoff(),
	     n_bigrams);

  if(recurs)
    for(INT i = 0; i < n_bigrams; i++)
      bigrams[i].print(f,l, noback);
}

inline VOID 
NGRAM_FINAL_OBJ::init(USHORT w, USHORT p, SHORT b) {
  set_word(w);
  set_prob(p);
  set_backoff(b);
}

USHORT 
NGRAM_FINAL_OBJ::read(GENFILE *infp, BOOL use_unk, USHORT unk_factor,
		      INT depth_left) {
  assert(depth_left == 0);
  if(infp->fread(&word, sizeof(USHORT), 1) < 1){
    perror("Reading ngram");
    //	panic("Ngram %d\n", ncreated);
  }
  if(infp->fread(&prob, sizeof(USHORT), 1) < 1){
    perror("Reading ngram");
    //	panic("Ngram %d\n", ncreated);
  }

  if(use_unk && word==Unknown_word_index)
    prob = MIN(65534, unk_factor + prob);

  //  if(Debug_flag.lm)
  //    printf("=NGRAM_FINAL_OBJ::read(): word = %d, prob =%d\n", word, prob);
      
  return(prob);
}

VOID 
NGRAM_FINAL_OBJ::write(GENFILE *outfp) {
  outfp->fwrite(&word, sizeof(USHORT), 1);
  outfp->fwrite(&prob, sizeof(USHORT), 1);
}

VOID
NGRAM_FINAL_OBJ::print(GENFILE *f, LEXICON *l, INT level, BOOL /*recurs*/) {
  for(INT i = 0; i < level; i++)
    f->fprintf("   ");

  f->fprintf("    %s(%d) %f(%d)  %f(%d)\n", l->vocabulary->get(word), word, 
	     logdecodelog(get_prob()), get_prob(),
	     logdecodelog(get_backoff()), get_backoff());
}

inline VOID 
NGRAM_INTER_OBJ::clear() {
  if(nextgrams != NULL){
    delete nextgrams;
    nextgrams = NULL;
  }
  //  n_nextgrams = 0;
}


inline NGRAM_FINAL_OBJ* 
NGRAM_INTER_OBJ::get_nextgram(INT w){
  UINT i, lower, upper;
  INT diff;
  
  if(nextgrams == NULL || get_n_nextgrams() == 0) return(NULL);

  lower = 0;
  upper = get_n_nextgrams();
  i = upper/2;
  while(upper > lower) {
    if((diff = w - nextgrams->get(i)->get_word()) == 0)
      return(nextgrams->get(i));
    else if(diff > 0)
      lower = i + 1;
    else // diff < 0
      upper = i;
    i = (lower + upper)/2;
  }
  return(NULL);
}

inline NGRAM_FINAL_OBJ* 
NGRAM_INTER_OBJ::get_nextgram_seq(INT w){
  static INT index = 0;
  NGRAM_FINAL_OBJ* res = NULL;
  INT diff;
  
  if(w == -1) { index = 0; return(NULL); }
  if(nextgrams == NULL || get_n_nextgrams() == 0) return(NULL);

  while(index < get_n_nextgrams()){
    if((diff = w - nextgrams->get(index)->get_word()) == 0) {
      res = nextgrams->get(index);
      index++; // static var
      break;
    }
    else if(diff > 0) {
      index++;
    }
    else { // diff < 0
      break;
    }
  }
  return(res);
}

VOID 
NGRAM_INTER_OBJ::load_nextgrams(INT n, NGRAM_ARRAY *ng) {
  if(n == 0) return;

  assert(nextgrams == NULL);
  //  n_nextgrams = n;
  if(ng->get_data()->ngram_mode() == FINAL){
    //    nextgrams = new NGRAM_FINAL_ARRAY(n_nextgrams);
    nextgrams = new NGRAM_FINAL_ARRAY(n);
  } else {
    //    nextgrams = new NGRAM_INTER_ARRAY(n_nextgrams);
    nextgrams = new NGRAM_INTER_ARRAY(n);
  }      
  nextgrams->load(ng, n);
}

USHORT 
NGRAM_INTER_OBJ::read(GENFILE *infp,  BOOL use_unk, USHORT unk_factor, 
		      INT depth_left){
  USHORT n_next;
  assert(nextgrams == NULL);
  assert(depth_left > 0);
  if(infp->fread(&word, sizeof(USHORT), 1) < 1){
    perror("Reading ngram");
    //	panic("Ngram %d\n", ncreated);
  }
  if(infp->fread(&prob, sizeof(USHORT), 1) < 1){
    perror("Reading ngram");
    //	panic("Ngram %d\n", ncreated);
  }
  if(infp->fread(&backoff, sizeof(SHORT), 1) < 1){
    perror("Reading ngram");
    //	panic("Ngram %d\n", ncreated);
  }
  if(infp->fread(&n_next, sizeof(USHORT), 1) < 1){
    perror("Reading ngram");
    //	panic("Ngram %d\n", ncreated);
  }

  if(use_unk && word==Unknown_word_index)
    prob = MIN(65534, unk_factor + prob);

  if(n_next > 0) {
    if(depth_left == 1){
      nextgrams = new NGRAM_FINAL_ARRAY(n_next);
    }
    else{
      nextgrams = new NGRAM_INTER_ARRAY(n_next);
    }

    USHORT p;
    depth_left--;
    for(INT i = 0; i < n_next; i++) 
      p = nextgrams->get(i)->read(infp, use_unk, unk_factor, depth_left);
  }

  //  if(Debug_flag.lm)
  //    printf("*NGRAM_INTER_OBJ::read(): word = %d, prob =%d, backoff = %d, nnext = %d\n", word, prob, backoff, n_nextgrams);

  return(prob);
}

VOID 
NGRAM_INTER_OBJ::write(GENFILE *outfp) {
  assert(get_n_nextgrams() <= (USHORT)0xffff);
  USHORT n_next = (USHORT)get_n_nextgrams();

  outfp->fwrite(&word, sizeof(USHORT), 1);
  outfp->fwrite(&prob, sizeof(USHORT), 1);
  outfp->fwrite(&backoff, sizeof(SHORT), 1);
  outfp->fwrite(&n_next, sizeof(USHORT), 1);

  for(INT i = 0; i < get_n_nextgrams(); i++)
    nextgrams->get(i)->write(outfp);
}

VOID
NGRAM_INTER_OBJ::print(GENFILE *f, LEXICON *l, INT level, BOOL recurs) {
  for(INT i = 0; i < level; i++)
    f->fprintf("   ");

  f->fprintf("    %s(%d) %f(%d)  %f(%d) (%d %d-grams)\n", 
	     l->vocabulary->get(word), word, 
	     logdecodelog(get_prob()), get_prob(),
	     logdecodelog(get_backoff()), get_backoff(),
	     get_n_nextgrams(), level+2);

  if(recurs)
    for(INT i = 0; i < get_n_nextgrams(); i++)
      nextgrams->get(i)->print(f,l,level+1,recurs);
}

INT
NGRAM::get_prob(LIST<INT>* history, INT nextword, INT hypid) {
  INT prob = PROB_NULL;

  if(use_unk) 
    context[0] = oov2unk(nextword);
  else
    context[0] = nextword;

  if(cached && hypid >= 0) {
    assert(cache_context[hypid] >= 0);
    prob = cache[hypid][context[0]];
  } 

  if(prob == PROB_NULL) {
    for(INT i = 1; i < order; i++)
      context[i] = Sentence_start_index;
    history->init_cursor();
    for(INT count = 1; count < order && history->next(); count++){
      // relies on left-right evaluation
      if(history->read() == Silence_index && history->next() == FALSE)
	break;
      if(use_unk)
	context[count] = oov2unk(history->read());
      else
	context[count] = history->read();
    }
    prob = get_prob_backoff(order);
    if(cached && hypid >= 0) 
      cache[hypid][context[0]] = prob;
  }
  return(prob);
}

INT
NGRAM::get_prob(HYP *h, INT nextword, INT hypid) {
  INT prob = PROB_NULL;
  HYP* p;

  if(use_unk) 
    context[0] = oov2unk(nextword);
  else
    context[0] = nextword;

  if(cached && hypid >= 0) {
    assert(cache_context[hypid] >= 0);
    prob = cache[hypid][context[0]];
  } 

  if(prob == PROB_NULL) {
    
    //    if(order == 2){
      //      if(use_unk) 
      //	context[1] = oov2unk(h->get_word());
      //      else
      //      context[1] = h->get_word();
    //      context[1] = h->get_lm_context()[0];
    //    } else 
    if(order > 1){
      memcpy(context+1, h->get_lm_context(), (order-1)*sizeof(INT));
    }

    //    if(use_unk) {
    //      for(INT i = 1; i < order; i++)
    //	context[i] = oov2unk(context[i]);
    //    }
	
    if(Debug_flag.lm)
      printf("prob = ");    
    prob = get_prob_backoff(order);
    if(Debug_flag.lm)
      printf("%d\n", prob);    
    if(cached && hypid >= 0) {
      cache[hypid][context[0]] = prob;
      if(Debug_flag.lm)
	printf("cache[%d][%d] = %d\n", hypid, nextword, prob);
    }
  } // prob == PROB_NULL
  else{
    if(Debug_flag.lm)
      printf("prob = cache[%d][%d] = %d\n", hypid, nextword, prob);  
  }
  h->lm_prob = prob;
  return(prob);
}

// Computing e.g. P(w_c0 | W_c2, W_c1) for trigram
//
// 1. P(w_n | w_n-(c-1) w_n-(c-2) ... w_n-1)
// 2. B(w_n-1, w_n-(c-1) w_n-(c-2) ... w_n-2)P(w_n |w_n-(c-2) ... w_n-1)
// 3. P(w_n |w_n-(c-2) ... w_n-1)
// 
INT
NGRAM::get_prob_backoff(INT clen) {
  assert(clen <= order);
  // root of the n-grams
  clen--;

  if(context[clen] >= get_nunigrams())
    return(get_prob_backoff(clen));

  NGRAM_FINAL_OBJ* ng = unigrams->get(context[clen]);

  if(clen == 0)
    return(ng->get_prob()*scale + offset);

  for(INT i = clen-1; i >= 1; i--) {
    ng = ng->get_nextgram(context[i]);
    if(ng == NULL)
      return(get_prob_backoff(clen));
  }
  NGRAM_FINAL_OBJ* fg = ng->get_nextgram(context[0]);
  if(fg == NULL){
    if(Debug_flag.lm)
      printf("%d + ", ng->get_backoff()*scale);    
    return(ng->get_backoff()*scale + get_prob_backoff(clen));  
  }  else {
    if(Debug_flag.lm)
      printf("%d", fg->get_prob()*scale + offset);
    return(fg->get_prob()*scale + offset);
  }
}

DOUBLE 
NGRAM::text_score(QUEUE<INT> *wq, GENFILE *outfp, LEXICON *lex, 
		  DOUBLE &noUNK_accum, INT &noov) {
  DOUBLE score = 0.0;
  LIST<INT> *cxt = new LIST<INT>();
  QUEUE<INT> *newq = new QUEUE<INT>();
  INT wd;
  DOUBLE oldp = noUNK_accum;

  if(Debug_flag.text_decode)
    fprintf(stderr, "In text_score()\n");

  // bbo default cxt is </s>
  // but for noway it is <s>
  for(INT j = 0; j < order-1; j++)
    cxt->push(Sentence_start_index);
    //    cxt->push(Sentence_end_index);

  while(wq->size > 0) {
    wd = wq->pop();
    if((wd != Sentence_start_index || wq->size==0) && wd != Paragraph_index && wd != Article_index) {
      INT p = get_prob(cxt, wd);
      //      score += p;
      score += logdecodelog(p);
      if(Verbose) {
	printf("log P(%s | ", lex->vocabulary->get(wd));
	cxt->init_cursor();
	for(INT count = 1; count < order && cxt->next(); count++)
	  printf("%s ", lex->vocabulary->get(cxt->read()));
	printf(") = %f (P = %f)\n", logdecodelog(p), logdecode(p));
      }
      if(wd != Unknown_word_index) {
	//	noUNK_accum += p;
	noUNK_accum += logdecodelog(p);
      } else {
	noov++;
      }
    }
    newq->push(wd);
    cxt->push(wd);
  }
  cxt->clear();
  outfp->fprintf(" (%.2f, %.2f including <UNK>)\n", noUNK_accum - oldp, 
	  score);
	  //	  logdecodelog(score)); 
  delete newq;
  delete cxt;
  return score;
}

BOOL 
NGRAM::hyp_similar(LIST<INT> *h1h, LIST<INT> *h2h){
  INT w1, w2;

  INT matched = 0;
  h1h->init_cursor();
  h2h->init_cursor();
  while(matched < order-1 && h1h->next() && h2h->next()) {
    w1 = h1h->read();
    w2 = h2h->read();
    if(w1 == Silence_index){
      if (h1h->next())
	w1 = h1h->read();
      else
	w1 = -1;
    }
    if(w2 == Silence_index){
      if(h2h->next())
	w2 = h2h->read();
      else
	w2 = -1;
    }
    if(w1 != w2) 
      return(FALSE);
    else
      matched++;
    if(w1 == -1)
      break;
  }
  if(matched == order-1 || (h1h->end_cursor() && h2h->end_cursor()))
    return(TRUE);
  return(FALSE);
}


BOOL
NGRAM::hyp_similar(HYP *h1, HYP *h2){
  for(INT i = 0; i < order-1; i++)
    if(h1->get_lm_context()[i] != h2->get_lm_context()[i])
      return(FALSE);
  return(TRUE);
}
  


BOOL 
NGRAM::hyp_similar(LIST<INT> *h1, LIST<INT> *h2, INT w){
  h2->push(w);
  BOOL res = hyp_similar(h1,h2);
  h2->pop();
  return(res);    
}

BOOL
NGRAM::hyp_similar(HYP *h1, HYP *h2, INT w){
  HYP* tmph;
  if(use_unk)
    tmph = new HYP(w, 0, 0, 0, h2, oov2unk(w));
  else
    tmph = new HYP(w, 0, 0, 0, h2);
  BOOL res = hyp_similar(h1, tmph);
  delete tmph;
  return(res);
}

INT
NGRAM::build(GENFILE *gramfp, LEXICON *lex, INT sc, INT os, BOOL cache_lm, 
	     INT nhy, BOOL build_vocab){
  LM_FORMAT fmt = NEWBIN;
  CHAR buf[256], *bptr;
  INT i, c;
  USHORT wd_idx;
  INT nuni;

  if(gramfp->fread(buf, sizeof(char), 4) != 4)
    panic("Failed to read LM header");

  scale = sc;
  offset = os;
  max_hyps = nhy;

  // check format and byte-swap status
  if(strcmp(buf, "TR1") == 0) {
    order = 3;
    fmt = OLDBIN;
  }
  else if(strcmp(buf, "TR2") == 0) {
    order = 3;
    fmt = NEWBIN;
  }
  else if(strcmp(buf, "2RT") == 0) {
    order = 3;
    fmt = NEWBIN;
    gramfp->flip_byteswap();
  } else if(strcmp(buf, "BI1") == 0) {
    order = 2;
    fmt = OLDBIN;
  }
  else if(strcmp(buf, "BI2") == 0) {
    panic("Bigram file created by noway [01].* must be used with -bigram (NOT -ngram)\n");
  }
  else if(strcmp(buf, "2IB") == 0) {
    panic("Bigram file created by noway [01].* must be used with -bigram (NOT -ngram)\n");
  } else if(strcmp(buf, "BI1") == 0) {
    panic("Bigram file created by noway [01].* must be used with -bigram (NOT -ngram)\n");
  } else if(strncmp(buf, "NG", 2) == 0 && buf[3] == 0) {
    fmt = NEWBIN;
    if(strpbrk(buf, "0123456789") == buf+2)
      order = atoi(buf+2);
    else
      panic("%s is unknown header to LM file\n", buf);      
  } else{
    // order not set yet
    fmt = ARPA;
  }

  if(Debug_flag.lm)
    printf("Building %d-gram LM\n", order);

  if(fmt != ARPA) {
    context = new INT[order];
    ngram_count = new INT[order];
  }
    
  // get nunigrams (and maybe nbigrams, ntrigrams)
  if(fmt == ARPA){
    if(Verbose){
      printf("ASCII ARPA format for ngram LM file\n");
      fflush(stdout);
    }
    // cheap hack to check if we have read any of the "/data/" from the file
    INT ch;
    if(strncmp(buf, arpa_hdr_chars, 4) == 0)
      ch = 4;
    else if(strncmp(&buf[1], arpa_hdr_chars, 3) == 0)
      ch = 3;
    else if(strncmp(&buf[2], arpa_hdr_chars, 2) == 0)
      ch = 2;
    else if(strncmp(&buf[3], arpa_hdr_chars, 1) == 0)
      ch = 1;
    else
      ch = 0;
    // still need to set order
    return(read_arpa_lm(gramfp, lex, ch, cache_lm,build_vocab));    
  } else { // fmt != ARPA
      assert(fmt != OLDBIN);
      for(INT level = 0; level < order; level++){
	if(gramfp->fread(&ngram_count[level], sizeof(INT), 1) != 1)
	  panic("LM failed to read n_%d-grams", level+1);
	if(Debug_flag.lm)
	  printf("%d %d-grams\n", ngram_count[level], level+1);
      }
      nuni = ngram_count[0];
  }

  // init various things
  if(order > 1)
    unigrams = new NGRAM_INTER_ARRAY(nuni);
  else
    unigrams = new NGRAM_FINAL_ARRAY(nuni);

  if(unk_factor == LM_NULL)
    unk_factor = (USHORT)MIN(50000,-logcode(1.0 / (FLOAT)nuni));


  if(cache_lm && max_hyps > 0 && order > 1)  
    create_cache();

  // read in LM lexicon
  // Note that it is okay for <UNK> to be in the LM even if 
  //  use_unk == FALSE - since the LM just scores word sequences
  //  (ie does not hypothesise following words) there is no problem
  //  with <UNK> being in the LM even if we do not use it.
  for(i = 0; i < get_nunigrams(); i++) {
    bptr = buf;
    while((c = gramfp->fgetc()) > 0){
      *bptr = (UCHAR)c;
      bptr++;
    }
    *bptr = '\0';
    gramfp->fread(&wd_idx, sizeof(USHORT), 1);
    if(build_vocab) {
      if(wd_idx != (USHORT)lex->vocabulary->size)
	panic("%s: wd_idx (%d) != vocab_size (%d)\n", buf, wd_idx, lex->vocabulary->size);
      lex->vocab_hash->insert(strdup(buf), wd_idx);
      lex->vocabulary->push(strdup(buf));

      if(!strcmp(buf, Sentence_start_symbol))
	Sentence_start_index = wd_idx;
      else if(!strcmp(buf, Sentence_end_symbol))
	Sentence_end_index = wd_idx;
      else if(!strcmp(buf, Unknown_word_symbol))
	Unknown_word_index = wd_idx;
      else if(!strcmp(buf, Paragraph_symbol))
	Paragraph_index = wd_idx;
      else if(!strcmp(buf, Article_symbol))
	Article_index = wd_idx;
      else if(!strcmp(buf, Silence_symbol))
	Silence_index = wd_idx;
      else ;
    } else {
      INT lexind = lex->vocab_hash->get(buf);
      assert(lexind == wd_idx);

      if(!strcmp(buf, Sentence_start_symbol))
	assert(Sentence_start_index == wd_idx);
      else if(!strcmp(buf, Sentence_end_symbol))
	assert(Sentence_end_index == wd_idx);
      else if(!strcmp(buf, Unknown_word_symbol))
	assert(Unknown_word_index == wd_idx);
      else if(!strcmp(buf, Paragraph_symbol))
	assert(Paragraph_index == wd_idx);
      else if(!strcmp(buf, Article_symbol))
	assert(Article_index == wd_idx);
      else if(!strcmp(buf, Silence_symbol))
	assert(Silence_index == wd_idx);
      else ;
    }
  }

  if(build_vocab){
    if(Sentence_start_index < 0)
      panic("Sentencestart symbol %s not in LM file", Sentence_start_symbol);

    if(Sentence_end_index < 0)
      Sentence_end_index = Sentence_start_index;
    //     panic("Sentence end symbol %s not in LM file", Sentence_end_symbol);

    if(lex->vocabulary->size > USHRT_MAX)
      panic("Vocabulary size too large: %d > %d", lex->vocabulary->size, 
	    USHRT_MAX);
  }
  assert(get_nunigrams() == lex->vocabulary->size);

  if(Verbose) {
    if(build_vocab)
      printf("Created vocabulary hash table (%d entries)\n", lex->vocab_hash->size());
    printf("Read %d unigrams\n", get_nunigrams());
    fflush(stdout);
  }

  if(fmt == OLDBIN) {
    if(Verbose){
      printf("Old binary format for trigram LM file\n");
      fflush(stdout);
    }
    return(read_old_bin_lm(gramfp, lex));
  }

  // read in unigrams and bigrams and trigrams
  if(Debug_flag.lm)
    printf("Reading %d unigrams\n", get_nunigrams());
  for(i = 0; i < get_nunigrams(); i++) {
    if(Debug_flag.lm)
      printf("Reading unigrams %d (at %ld)\n", i, gramfp->ftell());
    unigrams->get(i)->read(gramfp, use_unk, unk_factor, order-1);
    //    // Ugly hack that will get taken out when LMs conform to noway spec
    //    unigrams->get(i)->set_word(i);
    lex->unigrams->push(unigrams->get(i)->get_prob());
    //    lex->vocab_sorted->insert(new SCORE(i, unigrams[i].get_prob()));
  }

  if(gramfp->fgetc() != EOF)
    panic("Didn't reach end of file in LM (at %d)\n", gramfp->ftell());

  // check read expected numbers of unigrams/bigrams/trigrams
//  if(nunigrams != UNIGRAM_OBJ::ncreated)
//    panic("Expecting %d unigrams, read in %d\n", nunigrams,
//	  UNIGRAM_OBJ::ncreated);
//  if(nbigrams != BIGRAM_OBJ::ncreated)
//    panic("Expecting %d bigrams, read in %d\n", nbigrams,
//	  BIGRAM_OBJ::ncreated);
//  if(ntrigrams != TRIGRAM_OBJ::ncreated)
//    panic("Expecting %d trigrams, read in %d\n", ntrigrams,
//	  TRIGRAM_OBJ::ncreated);


  if(Verbose){
    printf("Built language model\n");
    info();
    fflush(stdout);
  }

  return(get_nunigrams());  
}

VOID
NGRAM::write_header(GENFILE *f, LEXICON *l){
  // write header
  CHAR hdr[4];
  sprintf(hdr, "NG%d", order);
  assert(strlen(hdr) == 3 && hdr[3] == 0);
  f->fwrite((void*)hdr, sizeof(CHAR), 4); // written as BIG_ENDIAN

  // write counts
  INT ng_count;
  for(INT level = 0; level < order; level++) {
    ng_count = ngram_count[level];
    f->fwrite(&ng_count, sizeof(INT), 1);
  }    
  
  // write out lexicon
  USHORT si;
  STRING wd;
  for(INT i = 0; i < get_nunigrams(); i++) {
    si = (USHORT) i;
    wd = l->vocabulary->get(i);
    f->fwrite((void*)wd, sizeof(CHAR), strlen(wd)+1);
    f->fwrite(&si, sizeof(USHORT), 1);
  }
}

VOID
NGRAM::write(GENFILE *f, LEXICON *l){
  write_header(f, l);
  for(INT j = 0; j < get_nunigrams(); j++)
    unigrams->get(j)->write(f);
}

VOID
NGRAM::print(GENFILE *f, LEXICON *l){
  for(INT i = 0; i < get_nunigrams(); i++)
    unigrams->get(i)->print(f, l, 0, TRUE);
}

inline VOID 
NGRAM::print_unigrams(GENFILE *f, LEXICON *l) {
  for(INT i = 0; i < get_nunigrams(); i++)
    unigrams->get(i)->print(f, l, 0, FALSE);
}


VOID 
NGRAM::info(GENFILE *outfp){
  fprintf(outfp->get_fp(), "%d-GRAM LANGUAGE MODEL:\n", order);
  for(INT i = 0; i < order; i++)
    fprintf(outfp->get_fp(), "%d %d-grams\n", ngram_count[i], i+1);
}

INT
NGRAM::read_arpa_lm(GENFILE *infp, LEXICON *lex, INT ch, BOOL cache_lm,
		    BOOL build_vocab){
  UINT i;
  USHORT si, sp;
  USHORT *wbuffer;
  SHORT sb=0;
  FLOAT prob, backoff;
  CHAR wd[256];
  CHAR buf[10001];
  INT nuni;

  // No data until we see "\data\" on a line
  if(ch > 0) {
    if(infp->eof())
      panic("LM file in unknown format (not BI1, BI2, TR1, TR2 or ARPA)");
    infp->fgets(buf, 10000);
    if(strncmp(buf, &arpa_hdr_chars[ch], strlen(&arpa_hdr_chars[ch])) != 0) {
      panic("LM file in unknown fmt (not BI1, BI2, TR1, TR2 or ARPA)");
    }
  } else { 
    do{
      if(infp->eof())

	panic("LM file in unknown format (not BI1, BI2, TR1, TR2 or ARPA)");
      infp->fgets(buf, 10000);
    } while(strncmp(buf, arpa_hdr_chars, strlen(arpa_hdr_chars)));
  }

  order = 0;
  ngram_count = new INT[16];
  INT count, level;
  while(1) {
    infp->fgets(buf, 10000);
    if(!strncmp(buf, "\\1-grams:", 9))
      break;
    i = 0;
    while(i < strlen(buf) && isspace(buf[i]))
      i++;
    CHAR *cptr = buf+i;
    if(i < strlen(buf)){ 
      if(sscanf(cptr, "ngram %d = %d\n", &level, &count) != 2)
	panic("Failed to read ARPA LM hdr:  %s\n (i = %d)", buf, i);
      order++;
      if(order > 16)
	panic("read_arpa_lm(): can't handle more than 15-grams!");
      if(order != level)
	panic("read_arpa_lm(): expecting ngram %d\n%s", order, cptr);
      ngram_count[order-1] = count;
      if(Verbose)
	printf("%d %d-grams\n", count, order);
    } // else line of whitespace      
  }
  assert(order > 0);
  context = new INT[order];
  wbuffer = new USHORT[order];
  nuni = ngram_count[0];
  if(order > 1)
    unigrams = new NGRAM_INTER_ARRAY(nuni);
  else
    unigrams = new NGRAM_FINAL_ARRAY(nuni);

  if(unk_factor == LM_NULL)
    unk_factor = (USHORT)MIN(50000,-logcode(1.0 / (FLOAT)nuni));

  BOOL contloop = FALSE;
  // loop over ngram order
  for(level = 1; level <= order; level++) {
    // buf contains "\\(level)-grams:"
    INT adj_count = 0;
    NGRAM_ARRAY* ngram_buf = NULL;
    NGRAM_FINAL_OBJ* parent = NULL;
    NGRAM_FINAL_OBJ* last_parent = NULL;
    INT ng_count = 0;
    if(level > 1){
      if(order > level)
	ngram_buf = new NGRAM_INTER_ARRAY(get_nunigrams());
      else
	ngram_buf = new NGRAM_FINAL_ARRAY(get_nunigrams());
    } // else write into unigrams

    // loop over level-grams
    for(i = 0, si = 0; i < (UINT)ngram_count[level-1]; i++){
      if(fscanf(infp->get_fp(), "%f", &prob) != 1)
	panic("read_arpa_lm(): failed to read %dth %d-gram prob", i, level);
      sp = (USHORT) (MAX(prob, -7.999) * -log_scale + 0.5);
      if(level == 1) {
	// unigrams are special as they are used to build the vocab
	if(fscanf(infp->get_fp(), "%s", wd) != 1)
	  panic("read_arpa_lm(): failed to read %dth unigram", i);
	if(!strcmp(wd, Unknown_word_symbol) && use_unk == FALSE) {
	  // ignore unknown words
	  contloop = TRUE;
	} else {
	  if(build_vocab){
	    lex->vocab_hash->insert(strdup(wd), (INT)si);
	    lex->vocabulary->push(strdup(wd));
	    wbuffer[0] = si;
	    si++;
	 
	    // special words...
	    if(!strcmp(wd, Unknown_word_symbol))
	      Unknown_word_index = wbuffer[0];
	    else if(!strcmp(wd, Sentence_start_symbol))
	      Sentence_start_index = wbuffer[0];
	    else if(!strcmp(wd, Sentence_end_symbol))
	      Sentence_end_index = wbuffer[0];
	    else if(!strcmp(wd, Paragraph_symbol))
	      Paragraph_index = wbuffer[0];
	    else if(!strcmp(wd, Article_symbol))
	      Article_index = wbuffer[0];
	    else if(!strcmp(wd, Silence_symbol))
	      Silence_index = wbuffer[0];
	    else ;
	  } else {
	  INT in = lex->vocab_hash->get(wd);
	  assert(in == (INT)si);
	  wbuffer[0] = si;
	  si++;
	  // Assume that special words are okay...
	  }
	} 
      } else { // level > 1
	for(INT o = 1; o <= level; o++){
	  if(fscanf(infp->get_fp(), "%s", wd) != 1)
	    panic("read_arpa_lm(): failed to read %dth %d-gram word", o+1, level+1); 
	  if(!strncmp(wd, Unknown_word_symbol, 5) && use_unk == FALSE) {
	    // ignore unknown words
	    contloop = TRUE;
	  } else {
	    wbuffer[o-1] = (USHORT)lex->vocab_hash->get(wd);
	    //	    if(use_unk == FALSE)
	    //	      uflag = TRUE;
	  }
	}
      }
      if(order > level) {
	if(fscanf(infp->get_fp(), "%f", &backoff) != 1)
	  panic("read_arpa_lm(): failed to read %dth %d-gram backoff", i, level+1);
	sb = (SHORT) (MAX(MIN(backoff,3.999), -3.999) * log_scale + 0.5);
      }

      if(contloop == TRUE){
	adj_count++;
	contloop = FALSE;
	continue;
      }

      // prob is sp, backoff is sb, words are in wbuffer
      if(wbuffer[level-1] == Unknown_word_index)
	sp = MIN(65534, unk_factor + sp);    
      if(level == 1) {
	unigrams->get(wbuffer[0])->set_word(wbuffer[0]);
	unigrams->get(wbuffer[0])->set_prob(sp);
	if(order > 1)
	  unigrams->get(wbuffer[0])->set_backoff(sb);
	if(build_vocab)
	  lex->unigrams->push(unigrams->get(wbuffer[0])->get_prob());
      } else {
	// get parent from wbuffer[0..level-2]
	parent = unigrams->get(wbuffer[0]);
	for(INT l = 1; l < level-1; l++)
	  parent = parent->get_nextgram(wbuffer[l]);
	if(ng_count > 0 && parent != last_parent) {
	  last_parent->load_nextgrams(ng_count,ngram_buf);
	  ng_count = 0;
	} 
	ngram_buf->get(ng_count)->init(wbuffer[level-1], sp, sb);
	ng_count++;
	last_parent = parent;
      }
    }
    if(ng_count > 0) {
      parent->load_nextgrams(ng_count,ngram_buf);
      ng_count = 0;
    } 
    if(level > 1)
      delete ngram_buf;
    ngram_count[level-1] -= adj_count;
    if(level == 1){
      // reset nunigrams
      //      nunigrams = ngram_count[0];
      unigrams->resize(ngram_count[0]);
      if(Sentence_start_index < 0)
	panic("Sentence start symbol %s not in LM file", Sentence_start_symbol);
      
      if(Sentence_end_index < 0)
	Sentence_end_index = Sentence_start_index;
	  //  	panic("Sentence end symbol %s not in LM file", Sentence_end_symbol);
      
      if(lex->vocabulary->size > USHRT_MAX)
	panic("Vocabulary size too large: %d > %d", lex->vocabulary->size, 
	      USHRT_MAX);
      assert(get_nunigrams() == lex->vocabulary->size);
      if(Verbose) 
	printf("Created vocabulary hash table (%d entries)\n", 
	       lex->vocab_hash->size());
    } 

    if(Verbose){
      printf("Read %d %d-grams (ignored %d)\n", 
	     ngram_count[level-1], level, adj_count);
      fflush(stdout);
    }
    if(level < order) {
      CHAR compstr[32];
      sprintf(compstr, "\\%d-grams:", level+1);
      do{
	infp->fgets(buf, 10000);
	if(infp->eof())
	  panic("read_arpa_lm(): premature eof (level = %d)", level);
      } while(strncmp(buf, compstr, 9));
    }
  }

  if(fscanf(infp->get_fp(), "%s", buf) != 1)
    panic("Failed to read \\end\\ in grammar file.");
  if(strcmp(buf, "\\end\\"))
    panic("Failed to read \\end\\ in grammar file");

  if(cache_lm && max_hyps > 0 && order > 1) 
    create_cache();

  return(get_nunigrams());
}


INT
NGRAM::read_old_bin_lm(GENFILE *, LEXICON *){
  panic("read_old_bin_lm(): not yet implemented");
  return(-1);
}

VOID
NGRAM::create_cache() {
  assert(max_hyps > 0 && order > 1);

  cache = new INT*[max_hyps];
  tmp_cache = new INT*[max_hyps];
  cache_context = new INT*[max_hyps];
  tmp_cache_context = new INT*[max_hyps];
  for(INT i = 0; i < max_hyps; i++){
    cache[i] = new INT[get_nunigrams()];
    tmp_cache[i] = new INT[get_nunigrams()];
    cache_context[i] = new INT[order-1];
    tmp_cache_context[i] = new INT[order-1];
  }
  backoff_cache = new INT[max_hyps];
  tmp_backoff_cache = new INT[max_hyps];
    

  caching = TRUE;

  clear_cache();
}

VOID
NGRAM::clear_cache() {
  assert(caching == TRUE);

  for(INT i = 0; i < max_hyps; i++) {
    backoff_cache[i] = 0;
    memset(cache[i], 0x7f, get_nunigrams()*sizeof(INT)); 
    for(INT k = 0; k < order-1; k++) {
      cache_context[i][k] = -1;
    }
  }
  cache_size = 0;
  cached = FALSE;
}

VOID
NGRAM::update_cache(LIST<HYP*>* hlist) {
  INT **iiptr;
  INT *iptr;
  INT indx;
  NGRAM_INTER_OBJ *ngo;
  NGRAM_FINAL_OBJ *fgo;
  HYP* h;
  HYP *p;

  // elements of hist are INT[order-1]

  assert(caching == TRUE);
  assert(hlist->size <= max_hyps);
    
  iiptr = tmp_cache;
  tmp_cache = cache;
  cache = iiptr;
 
  iiptr = tmp_cache_context;
  tmp_cache_context = cache_context;
  cache_context = iiptr;
 
  iptr = tmp_backoff_cache;
  tmp_backoff_cache = backoff_cache;
  backoff_cache = iptr;

  tmp_cache_size = cache_size;
  cache_size = 0;

  if(Debug_flag.trigram)
    printf("\nBuilding Cache\n");
  hlist->init_cursor();
  while(hlist->next()) {
    if(Debug_flag.trigram)
      hlist->read()->info();

    h = hlist->read();

    if((indx = in_tmp_cache(h)) == -1) {
      // Compute a row of trigram probs
      memset(cache[cache_size], 0x7f, sizeof(INT)*get_nunigrams());
      if((ngo = get_context_obj(h)) == NULL) {
	backoff_cache[cache_size] = 0;
	if(Debug_flag.trigram)
	  printf("backoff_cache[%d] = %d\n", cache_size, backoff_cache[cache_size]);
      } else {
	backoff_cache[cache_size] = ngo->get_backoff()*scale;
	if(Debug_flag.trigram)
	  printf("backoff_cache[%d] = %d\n", cache_size, backoff_cache[cache_size]);
	for(INT j = 0; j < ngo->get_n_nextgrams(); j++){
	  fgo = ngo->get_nextgram_by_index(j);
	  cache[cache_size][fgo->get_word()] = fgo->get_prob()*scale + offset;
	  if(Debug_flag.trigram)
	    printf("cache[%d][%d] = %d\n", cache_size, fgo->get_word(),
		   cache[cache_size][fgo->get_word()]);
	}
      }
    } else {
      // Copy probs and backoff from previous cache
      if(Debug_flag.trigram)
	  printf("copy prev_cache[%d] to cache[%d]\n", indx, cache_size);
      iptr = cache[cache_size];
      cache[cache_size] = tmp_cache[indx];
      tmp_cache[indx] = iptr;
      backoff_cache[cache_size] = tmp_backoff_cache[indx];
    }

    // Finally update cache_context
    //    if(order == 2){
    //      cache_context[cache_size][0] = h->get_word();
    //      if(Debug_flag.trigram)
    //	printf("cache_context[%d][0] = %d\n", cache_size, 
    //	       cache_context[cache_size][0]);
    //    } else if(order == 3){
    //      cache_context[cache_size][0] = h->get_word();
    //      cache_context[cache_size][1] = h->get_prev_word();
    //      if(Debug_flag.trigram)
    //	printf("cache_context[%d][1] = %d\n", cache_size, 
    //	       cache_context[cache_size][1]);
    //    } else 

    if(order > 2){
      //      cache_context[cache_size][0] = h->get_word();
      memcpy(cache_context[cache_size], h->get_lm_context(), (order-1)*sizeof(INT));
    }

    if(Debug_flag.cache){
      fprintf(stderr, "cache_context[h=%d] = ", cache_size);
      for(INT k = order-2; k >= 0; k--) 
	fprintf(stderr, "%d ", cache_context[cache_size][k]);
      fprintf(stderr, "\n");
    }
    cache_size++;
  }
  cached = TRUE;
}


// returns -1 if context not in tmp_cache
INT 
NGRAM::in_tmp_cache(HYP *h){
  BOOL match;
  HYP *p;

  for(INT i = 0; i < tmp_cache_size; i++){
    match = TRUE;
    if(order > 1){
      //if (tmp_cache_context[i][0] != h->get_word()){
      //	match = FALSE;
      //	continue; // next i
      //      }
      //      if(order > 2){
      for(INT j = 0; j < order-1; j++){
	if (tmp_cache_context[i][j] != h->get_lm_context()[j]){
	  match = FALSE;
	  break; // out of loop j
	}
      }
    }
    //    }
    if(match) return(i);
  } // for i
  // No match
  return(-1);
}


NGRAM_INTER_OBJ* 
NGRAM::get_context_obj(HYP* h){
  HYP* p;
  NGRAM_INTER_OBJ* ng = NULL;

  if(order > 1){
    ng =  (NGRAM_INTER_OBJ*) unigrams->get(h->get_lm_context()[order-2]);
    for(INT i = order-3; i >= 0; i++){
      ng = (NGRAM_INTER_OBJ*) ng->get_nextgram(h->get_lm_context()[i]);
      if(ng == NULL) return NULL;
    }
  }
      
  return(ng);
}

INT
TRIGRAM::build(GENFILE *gramfp, LEXICON *lex, INT sc, INT os, BOOL cache_lm, 
	       INT nhy, BOOL build_vocab){
  LM_FORMAT fmt = NEWBIN;
  CHAR buf[256], *bptr;
  INT i, c;
  USHORT wd_idx;

  if(gramfp->fread(buf, sizeof(char), 4) != 4)
    panic("Failed to read LM header");

  scale = sc;
  offset = os;
  max_hyps = nhy;

  // check format and byte-swap status
  if(strcmp(buf, "TR1") == 0) {
    order = 3;
    fmt = OLDBIN;
  }
  else if(strcmp(buf, "TR2") == 0) {
    order = 3;
    fmt = NEWBIN;
  }
  else if(strcmp(buf, "2RT") == 0) {
    order = 3;
    fmt = NEWBIN;
    gramfp->flip_byteswap();
  } else if(strcmp(buf, "BI1") == 0) {
    order = 2;
    fmt = OLDBIN;
  }
  else if(strcmp(buf, "NG3") == 0) {
    fmt = NEWBIN;
    order = 3;
  }
  else if (strcmp(buf, "BI1") == 0 || strcmp(buf, "BI2") == 0 
	   || strcmp(buf, "1IB") == 0 || strcmp(buf, "2IB") == 0) {
    panic("Bigram LM (%s) but -trigram specified!\n", buf);
  } 
  else if (strncmp(buf, "NG", 2) == 0 && buf[3] == 0 && buf[2] != '3') {
    panic("NGRAM LM (order %c) but -trigram specified!\n", buf[2]);
  }
  else{
    // order not set yet
    fmt = ARPA;
  }

  if(Debug_flag.lm)
    printf("Building %d-gram LM\n", order);

  if(fmt != ARPA) {
    context = new INT[order];
    ngram_count = new INT[order];
  }
    
  // get nunigrams (and maybe nbigrams, ntrigrams)
  if(fmt == ARPA){
    panic("Use -ngram for ARPA format LMs\n");
  } else { // fmt != ARPA
      assert(fmt != OLDBIN);
      for(INT level = 0; level < order; level++){
	if(gramfp->fread(&ngram_count[level], sizeof(INT), 1) != 1)
	  panic("LM failed to read n_%d-grams", level+1);
	if(Debug_flag.lm)
	  printf("%d %d-grams\n", ngram_count[level], level+1);
      }
      set_nunigrams(ngram_count[0]);
  }


  // init various things
  tg_unigrams = new UNIGRAM_OBJ[get_nunigrams()];

  if(unk_factor == LM_NULL)
    unk_factor = (USHORT)MIN(50000,-logcode(1.0 / (FLOAT)get_nunigrams()));


  if(cache_lm && max_hyps > 0)
    create_cache();

  // read in LM lexicon
  // Note that it is okay for <UNK> to be in the LM even if 
  //  use_unk == FALSE - since the LM just scores word sequences
  //  (ie does not hypothesise following words) there is no problem
  //  with <UNK> being in the LM even if we do not use it.
  for(i = 0; i < get_nunigrams(); i++) {
    bptr = buf;
    while((c = gramfp->fgetc()) > 0){
      *bptr = (UCHAR)c;
      bptr++;
    }
    *bptr = '\0';
    gramfp->fread(&wd_idx, sizeof(USHORT), 1);
    if(wd_idx != (USHORT)lex->vocabulary->size)
      panic("%s: wd_idx (%d) != vocab_size (%d)\n", buf, wd_idx, lex->vocabulary->size);
    lex->vocab_hash->insert(strdup(buf), wd_idx);
    lex->vocabulary->push(strdup(buf));

    if(!strcmp(buf, Sentence_start_symbol))
      Sentence_start_index = wd_idx;
    else if(!strcmp(buf, Sentence_end_symbol))
      Sentence_end_index = wd_idx;
    else if(!strcmp(buf, Unknown_word_symbol))
      Unknown_word_index = wd_idx;
    else if(!strcmp(buf, Paragraph_symbol))
      Paragraph_index = wd_idx;
    else if(!strcmp(buf, Article_symbol))
      Article_index = wd_idx;
    else if(!strcmp(buf, Silence_symbol))
      Silence_index = wd_idx;
    else ;
  }

  if(Sentence_start_index < 0)
     panic("Sentence start symbol %s not in LM file", Sentence_start_symbol);

  if(Sentence_end_index < 0)
    Sentence_end_index = Sentence_start_index;
    //     panic("Sentence end symbol %s not in LM file", Sentence_end_symbol);

  if(lex->vocabulary->size > USHRT_MAX)
    panic("Vocabulary size too large: %d > %d", lex->vocabulary->size, 
	  USHRT_MAX);
  assert(get_nunigrams() == lex->vocabulary->size);

  if(Verbose) {
    printf("Created vocabulary hash table (%d entries)\n", lex->vocab_hash->size());
    printf("Read %d unigrams\n", get_nunigrams());
    fflush(stdout);
  }

  if(fmt == OLDBIN) {
    panic("Cannot handle old binary LM fmt!\n");
  }

  // read in unigrams and bigrams and trigrams
  if(Debug_flag.lm)
    printf("Reading %d unigrams\n", get_nunigrams());
  for(i = 0; i < get_nunigrams(); i++) {
    tg_unigrams[i].read(gramfp, use_unk, unk_factor);
    lex->unigrams->push(tg_unigrams[i].get_prob());
  }

  if(gramfp->fgetc() != EOF)
    panic("Didn't reach end of file in LM (at %d)\n", gramfp->ftell());

  if(Verbose){
    printf("Built language model\n");
    info();
    fflush(stdout);
  }

  return(get_nunigrams());  
}



VOID
TRIGRAM::update_cache(LIST<HYP*>* hlist) {
  INT **iiptr;
  INT *iptr;
  INT indx;
  BIGRAM_OBJ* bgo;
  TRIGRAM_OBJ* tgo;
  HYP* h;

  // elements of hist are INT[order-1]

  assert(caching == TRUE);
  assert(hlist->size <= max_hyps);
    
  iiptr = tmp_cache;
  tmp_cache = cache;
  cache = iiptr;
 
  iiptr = tmp_cache_context;
  tmp_cache_context = cache_context;
  cache_context = iiptr;
 
  iptr = tmp_backoff_cache;
  tmp_backoff_cache = backoff_cache;
  backoff_cache = iptr;

  tmp_cache_size = cache_size;
  cache_size = 0;

  if(Debug_flag.trigram)
    printf("\nBuilding Cache\n");
  hlist->init_cursor();
  while(hlist->next()) {
    if(Debug_flag.trigram)
      hlist->read()->info();

    h = hlist->read();

    if((indx = in_tmp_cache(h)) == -1) {
      // Compute a row of trigram probs
      memset(cache[cache_size], 0x7f, sizeof(INT)*get_nunigrams());
      if((bgo = get_bigram_obj(h)) == NULL) {
	backoff_cache[cache_size] = 0;
	if(Debug_flag.trigram)
	  printf("backoff_cache[%d] = %d\n", cache_size, backoff_cache[cache_size]);
      } else {
	backoff_cache[cache_size] = bgo->get_backoff()*scale;
	if(Debug_flag.trigram)
	  printf("backoff_cache[%d] = %d\n", cache_size, backoff_cache[cache_size]);
	for(INT j = 0; j < bgo->get_n_trigrams(); j++){
	  tgo = bgo->get_trigram_by_index(j);
	  cache[cache_size][tgo->get_word()] = tgo->get_prob()*scale + offset;
	  if(Debug_flag.trigram)
	    printf("cache[%d][%d] = %d\n", cache_size, tgo->get_word(),
		   cache[cache_size][tgo->get_word()]);
	}
      }
    } else {
      // Copy probs and backoff from previous cache
      if(Debug_flag.trigram)
	  printf("copy prev_cache[%d] to cache[%d]\n", indx, cache_size);
      iptr = cache[cache_size];
      cache[cache_size] = tmp_cache[indx];
      tmp_cache[indx] = iptr;
      backoff_cache[cache_size] = tmp_backoff_cache[indx];
    }

    // Finally update cache_context
    cache_context[cache_size][0] = h->get_word();
    if(Debug_flag.trigram)
	printf("cache_context[%d][0] = %d\n", cache_size, 
	       cache_context[cache_size][0]);
    cache_context[cache_size][1] = h->get_prev_word();
    if(Debug_flag.trigram)
      printf("cache_context[%d][1] = %d\n", cache_size, 
	     cache_context[cache_size][1]);

    if(Debug_flag.cache){
      fprintf(stderr, "cache_context[h=%d] = ", cache_size);
      for(INT k = order-2; k >= 0; k--) 
	fprintf(stderr, "%d ", cache_context[cache_size][k]);
      fprintf(stderr, "\n");
    }
    cache_size++;
  }
  cached = TRUE;
}



BOOL 
TRIGRAM::hyp_similar(LIST<INT> *h1h, LIST<INT> *h2h) {
  INT w1, w2;

  h1h->init_cursor();
  h2h->init_cursor();
  INT matched = 0;
  while(matched < 2 && h1h->next() && h2h->next()) {
    w1 = h1h->read();
    w2 = h2h->read();
    if(w1 == Silence_index){
      if (h1h->next())
	w1 = h1h->read();
      else
	w1 = -1;
    }
    if(w2 == Silence_index){
      if(h2h->next())
	w2 = h2h->read();
      else
	w2 = -1;
    }
    if(w1 != w2) 
      return(FALSE);
    else
      matched++;
    if(w1 == -1)
      break;
  }
  if(matched == 2 || (h1h->end_cursor() && h2h->end_cursor()))
    return(TRUE);
  return(FALSE);
}

BOOL 
TRIGRAM::hyp_similar(LIST<INT> *h1h, LIST<INT> *h2h, INT w) {
  INT w1, w2;

  h1h->init_cursor();
  h2h->init_cursor();

  // First word
  if(h1h->next())
    w1 = h1h->read();
  else
    w1 = -1;

  if(w1==Silence_index) {
    if(h1h->next())
      w1 = h1h->read();
    else
      w1 = -1;
  }
      
  if(w1 != w)
    return(FALSE);

  if(h1h->next())
    w1 = h1h->read();
  else
    w1 = -1;

  if(w1==Silence_index) {
    if(h1h->next())
      w1 = h1h->read();
    else
      w1 = -1;
  }
      
  if(h2h->next())
    w2 = h2h->read();
  else
    w2 = -1;

  if(w2==Silence_index) {
    if(h2h->next())
      w2 = h2h->read();
    else
      w2 = -1;
  }
      
  return(w1 == w2);
}


INT
TRIGRAM::get_prob(LIST<INT>* history, INT w3, INT hypid) {
  INT prob;
  INT w1, w2;

  if(use_unk) 
    w3 = oov2unk(w3);

  if(cached && hypid >= 0) {
    prob = cache[hypid][w3];
    if(prob == PROB_NULL) {
      // we know that we have to backoff to get p(w | hist)
      // interword <SIL> is never top
      if (history->size == 0) {
	w2 = Sentence_start_index;
      } else {
	w2 = history->top();
	if(use_unk) w2 = oov2unk(w2);
      }
      prob = get_prob(w2, w3) + backoff_cache[hypid];
      cache[hypid][w3] = prob; 
    }
  } else {
    // interword <SIL> is never top
    if (history->size == 0) {
      w1 = Sentence_start_index;
      w2 = Sentence_start_index;
    } else {
      w2 = history->top();
      if(use_unk) w2 = oov2unk(w2);
      if(history->size > 1) {
	w1 = history->nth(2);
	if(w1 == Silence_index) {
	  if(history->size > 2)
	    w1 = history->nth(3);
	  else
	    w1 = Sentence_start_index;
	}
	if(use_unk) w1 = oov2unk(w1);
      }	else {
	w1 = Sentence_start_index;
      }
    }
    prob = get_prob(w1, w2, w3);
  }
  return(prob);
}

INT
TRIGRAM::get_prob(HYP* h,  INT w3, INT hypid) {
  INT w1, w2;
  INT prob; 
  if(use_unk) 
    w3 = oov2unk(w3);

  if(cached && hypid >= 0) {
    prob = cache[hypid][w3];
    if(prob == PROB_NULL) {
      //      if(use_unk)
      //	w2 = oov2unk(h->get_word());
      //      else
      //	w2 = h->get_word();
      w2 = h->get_lm_context()[0];
      prob = get_prob(w2, w3) + backoff_cache[hypid];
      if(Debug_flag.trigram)
	printf("prob = get_prob(%d, %d) + backoff_cache[%d] = %d+%d = %d\n", 
	       w2, w3, hypid, get_prob(w2, w3), 
	       backoff_cache[hypid], prob);      
      cache[hypid][w3] = prob; 
    }else {
      if(Debug_flag.trigram)
	printf("prob = cache[%d][%d] = %d\n", hypid, w3, prob);
    }
  } else {
    //    w1 = h->get_prev_word();
    //    w2 = h->get_word();
    //    if(use_unk){
    //      w1 = oov2unk(w1);
    //      w2 = oov2unk(w2);
    //    }
    //    prob = get_prob(w1, w2, w3);
    prob = get_prob(h->get_lm_context()[1], (h->get_lm_context()[0], w3));
  }
  h->lm_prob = prob;
  return(prob);
}

INT
TRIGRAM::get_prob(INT w1, INT w2, INT w3){
  BIGRAM_OBJ *bi;
  TRIGRAM_OBJ *tri = NULL;

  if((bi = tg_unigrams[w1].get_bigram(w2)) != NULL) {
    if(((tri = bi->get_trigram(w3)) != NULL)) {
      return(tri->get_prob()*scale + offset);
    } else {
      // offset in get_prob(w2,w3)
      return(bi->get_backoff()*scale + get_prob(w2, w3));
    }
  }
  else{
    // offset in get_prob(w2,w3)
    return(get_prob(w2, w3));
  }
  return(0);
}

INT
TRIGRAM::get_prob(INT w1, INT w2) {
    BIGRAM_OBJ *bi;

    if((bi = tg_unigrams[w1].get_bigram(w2)) != NULL) {
      return(bi->get_prob()*scale + offset);
    } else {
      return((tg_unigrams[w1].get_backoff() + 
	      tg_unigrams[w2].get_prob())*scale + offset);
    }
    return(0);
}

VOID
TRIGRAM::write(GENFILE *f, LEXICON *l){
  write_header(f, l);
  for(INT j = 0; j < get_nunigrams(); j++)
    tg_unigrams[j].write(f);
}



INT
BIGRAM::build(GENFILE *gramfp, LEXICON *lex, INT sc, INT os, BOOL cache_lm, 
	      INT nhy, BOOL build_vocab){
  LM_FORMAT fmt = NEWBIN;
  CHAR buf[256], *bptr;
  INT i, c;
  USHORT wd_idx;

  if(gramfp->fread(buf, sizeof(char), 4) != 4)
    panic("Failed to read LM header");

  scale = sc;
  offset = os;
  max_hyps = nhy;

  // check format and byte-swap status
  if(strcmp(buf, "BI1") == 0) {
    order = 2;
    fmt = OLDBIN;
  }
  else if(strcmp(buf, "BI2") == 0) {
    order = 2;
    fmt = NEWBIN;
  }
  else if(strcmp(buf, "2IB") == 0) {
    order = 2;
    fmt = NEWBIN;
    gramfp->flip_byteswap();
  } 
  else if(strcmp(buf, "NG2") == 0){
    order = 2;
    fmt = NEWBIN;
  } 
  else if (strcmp(buf, "TR1") == 0 || strcmp(buf, "TR2") == 0 
	   || strcmp(buf, "1RT") == 0 || strcmp(buf, "2RT") == 0) {
    panic("Trigram LM (%s) but -bigram specified!\n", buf);
  } 
  else if (strncmp(buf, "NG", 2) == 0 && buf[3] == 0 && buf[2] != '2') {
    panic("NGRAM LM (order %c) but -bigram specified!\n", buf[2]);
  }
  else{
    // order not set yet
    fmt = ARPA;
  }

  if(Debug_flag.lm)
    printf("Building %d-gram LM\n", order);

  if(fmt != ARPA) {
    context = new INT[order];
    ngram_count = new INT[order];
  }
    
  // get nunigrams (and maybe nbigrams, ntrigrams)
  if(fmt == ARPA){
    panic("Use -ngram (not -bigram) for ARPA format LMs\n");
  } else { // fmt != ARPA
      assert(fmt != OLDBIN);
      for(INT level = 0; level < order; level++){
	if(gramfp->fread(&ngram_count[level], sizeof(INT), 1) != 1)
	  panic("LM failed to read n_%d-grams", level+1);
	if(Debug_flag.lm)
	  printf("%d %d-grams\n", ngram_count[level], level+1);
      }
      set_nunigrams(ngram_count[0]);
  }


  // init various things
  bg_unigrams = new UNIGRAM_OBJ[get_nunigrams()];

  if(unk_factor == LM_NULL)
    unk_factor = (USHORT)MIN(50000,-logcode(1.0 / (FLOAT)get_nunigrams()));


  if(cache_lm && max_hyps > 0)
    create_cache();

  // read in LM lexicon
  // Note that it is okay for <UNK> to be in the LM even if 
  //  use_unk == FALSE - since the LM just scores word sequences
  //  (ie does not hypothesise following words) there is no problem
  //  with <UNK> being in the LM even if we do not use it.
  for(i = 0; i < get_nunigrams(); i++) {
    bptr = buf;
    while((c = gramfp->fgetc()) > 0){
      *bptr = (UCHAR)c;
      bptr++;
    }
    *bptr = '\0';
    gramfp->fread(&wd_idx, sizeof(USHORT), 1);
    if(wd_idx != (USHORT)lex->vocabulary->size)
      panic("%s: wd_idx (%d) != vocab_size (%d)\n", buf, wd_idx, lex->vocabulary->size);
    lex->vocab_hash->insert(strdup(buf), wd_idx);
    lex->vocabulary->push(strdup(buf));

    if(!strcmp(buf, Sentence_start_symbol))
      Sentence_start_index = wd_idx;
    else if(!strcmp(buf, Sentence_end_symbol))
      Sentence_end_index = wd_idx;
    else if(!strcmp(buf, Unknown_word_symbol))
      Unknown_word_index = wd_idx;
    else if(!strcmp(buf, Paragraph_symbol))
      Paragraph_index = wd_idx;
    else if(!strcmp(buf, Article_symbol))
      Article_index = wd_idx;
    else if(!strcmp(buf, Silence_symbol))
      Silence_index = wd_idx;
    else ;
  }

  if(Sentence_start_index < 0)
     panic("Sentence start symbol %s not in LM file", Sentence_start_symbol);

  if(Sentence_end_index < 0)
    Sentence_end_index = Sentence_start_index;
    //     panic("Sentence end symbol %s not in LM file", Sentence_end_symbol);

  if(lex->vocabulary->size > USHRT_MAX)
    panic("Vocabulary size too large: %d > %d", lex->vocabulary->size, 
	  USHRT_MAX);
  assert(get_nunigrams() == lex->vocabulary->size);

  if(Verbose) {
    printf("Created vocabulary hash table (%d entries)\n", lex->vocab_hash->size());
    printf("Read %d unigrams\n", get_nunigrams());
    fflush(stdout);
  }

  if(fmt == OLDBIN) {
    panic("Cannot handle old binary LM fmt!\n");
  }

  // read in unigrams and bigrams and trigrams
  if(Debug_flag.lm)
    printf("Reading %d unigrams\n", get_nunigrams());
  for(i = 0; i < get_nunigrams(); i++) {
    bg_unigrams[i].read(gramfp, use_unk, unk_factor);
    lex->unigrams->push(bg_unigrams[i].get_prob());
  }

  if(gramfp->fgetc() != EOF)
    panic("Didn't reach end of file in LM (at %d)\n", gramfp->ftell());

  if(Verbose){
    printf("Built language model\n");
    info();
    fflush(stdout);
  }

  return(get_nunigrams());  
}



VOID
BIGRAM::update_cache(LIST<HYP*>* hlist) {
  INT **iiptr;
  INT *iptr;
  INT indx;
  UNIGRAM_OBJ* ugo;
  BIGRAM_OBJ* bgo;
  HYP* h;

  // elements of hist are INT[order-1]

  assert(caching == TRUE);
  assert(hlist->size <= max_hyps);
    
  iiptr = tmp_cache;
  tmp_cache = cache;
  cache = iiptr;
 
  iiptr = tmp_cache_context;
  tmp_cache_context = cache_context;
  cache_context = iiptr;
 
  iptr = tmp_backoff_cache;
  tmp_backoff_cache = backoff_cache;
  backoff_cache = iptr;

  tmp_cache_size = cache_size;
  cache_size = 0;

  hlist->init_cursor();
  while(hlist->next()) {
    h = hlist->read();

    if((indx = in_tmp_cache(h)) == -1) {
      // Compute a row of trigram probs
      memset(cache[cache_size], 0x7f, sizeof(INT)*get_nunigrams());
      ugo = &(bg_unigrams[h->get_word()]);
      backoff_cache[cache_size] = ugo->get_backoff();
      for(INT j = 0; j < ugo->get_n_bigrams(); j++){
	bgo = ugo->get_bigram_by_index(j);
	cache[cache_size][bgo->get_word()] = bgo->get_prob()*scale + offset;
      }
    } else {
      // Copy probs and backoff from previous cache
      iptr = cache[cache_size];
      cache[cache_size] = tmp_cache[indx];
      tmp_cache[indx] = iptr;
      backoff_cache[cache_size] = tmp_backoff_cache[indx];
    }

    // Finally update cache_context
    cache_context[cache_size][0] = h->get_word();
    cache_size++;
  }
  cached = TRUE;
}


INT
BIGRAM::get_prob(LIST<INT>* history, INT w2, INT hypid) {
  INT w1;
  INT prob; 

  if(use_unk) 
    w2 = oov2unk(w2);

  if(cached && hypid >= 0) {
    prob = cache[hypid][w2];
    if(prob == PROB_NULL) {
      if(use_unk)
	w1 = oov2unk(history->top());
      else
	w1 = history->top();
      prob = (bg_unigrams[w2].get_prob() + backoff_cache[hypid])*scale + offset;
      cache[hypid][w2] = prob; 
    }
  } else {
    w1 = history->top();
    if(use_unk)
      w1 = oov2unk(w1);
    prob = get_prob(w1, w2);
  }
  return(prob);
}


INT
BIGRAM::get_prob(HYP* h,  INT w2, INT hypid) {
  INT w1;
  INT prob; 
  if(use_unk) 
    w2 = oov2unk(w2);

  if(cached && hypid >= 0) {
    prob = cache[hypid][w2];
    if(prob == PROB_NULL) {
      //      if(use_unk)
      //	w1 = oov2unk(h->get_word());
      //      else
      //	w1 = h->get_word();
      w1 = h->get_lm_context()[0];
      prob = (bg_unigrams[w2].get_prob() + backoff_cache[hypid])*scale + offset;
      cache[hypid][w2] = prob; 
    }
  } else {
    //    w1 = h->get_word();
    //    if(use_unk)
    //      w1 = oov2unk(w1);
    prob = get_prob(h->get_lm_context()[0], w2);
  }
  h->lm_prob = prob;
  return(prob);
}

INT
BIGRAM::get_prob(INT w1, INT w2) {
    BIGRAM_OBJ *bi;

    if((bi = bg_unigrams[w1].get_bigram(w2)) != NULL) {
      return(bi->get_prob()*scale + offset);
    } else {
      return((bg_unigrams[w1].get_backoff() + 
	      bg_unigrams[w2].get_prob())*scale + offset);
    }
    return(0);
}


VOID
BIGRAM::write(GENFILE *f, LEXICON *l){
  write_header(f, l);
  for(INT j = 0; j < get_nunigrams(); j++)
    bg_unigrams[j].write(f, TRUE);
}

