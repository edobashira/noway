// -*- Mode: C++;  -*-
// File: NW_mixture_ngram.cc
// Author: Steve Renals (sjr@dcs.shef.ac.uk)
// Copyright (C) Department of Computer Science, Sheffield University, 1997-98
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

#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_debug.h"
#include "NW_param.h"

#include "NW_hypothesis.h"
#include "NW_mixture_ngram.h"
#include "NW_lexicon.h"

// File format:
//
// 
// <ngram-file> <floating-point-weight>

INT
MIXTURE_NGRAM::build(GENFILE *f, LEXICON *l, INT sc, INT os, 
			  BOOL cache_lm, INT nhy, BOOL build_vocab) {
  CHAR buf[BUFSIZ+1];
  CHAR fname[BUFSIZ+1];
  GENFILE* ng_fp;
  NGRAM* ng;
  INT nu;
  
  
  assert(cache_lm == FALSE);

  scale = sc;
  offset = os;
  max_hyps = nhy;

  while(f->fgets(buf, BUFSIZ) != NULL) {
    if(sscanf(buf, "%s", fname) != 1)
      continue;
    ng_fp = new GENFILE(fname, "r");
    if(Verbose)
      printf("%s\n", fname);
    // no scale, offset or caching for component LMs
    ng = new NGRAM();
    nu = ng->build(ng_fp, l, 1, 0, FALSE, max_hyps, build_vocab);
    build_vocab = FALSE;
    if(get_nunigrams() == -1)
      set_nunigrams(nu);
    else
      assert(nu == get_nunigrams());
    //    ng_probs->push(0);
    ng_lst->push(ng);
    order = MAX(order, ng_lst->top()->get_order());
    delete ng_fp;
  }

  uniform_wt = logcode(1.0/(FLOAT)(ng_lst->size));

  if(cache_lm && max_hyps > 0)
    create_cache();

  return(get_nunigrams());
}


INT
LSA_MIXTURE_NGRAM::build(GENFILE *f, LEXICON *l, INT sc, INT os, 
			  BOOL cache_lm, INT nhy, BOOL build_vocab) {
  return 0;
}

VOID 
MIXTURE_NGRAM::update_cache(LIST<HYP*>* hlist){
}


INT
MIXTURE_NGRAM::get_prob(HYP* h, INT nextword, INT hypid){
  get_component_probs(h, nextword, hypid); 
  panic("Not yet implemented!\n");
  return h->wt_sum_probs();
}

INT 
BLIND_MIXTURE_NGRAM::get_prob(HYP* h, INT nextword, INT hypid){
  get_component_probs(h, nextword, hypid); // sets h->ng_probs_i to tg_i(n|n-1,n-2)
  // h->update_ng_wts();
  return h->wt_sum_probs();  //sets h->lm_prob  and wts ng_probs
}
    
INT 
LSA_MIXTURE_NGRAM::get_prob(HYP* h, INT nextword, INT hypid){
  return 0;
}

VOID
MIXTURE_NGRAM::get_component_probs(HYP* h, INT nextword, INT hypid){
  ng_lst->init_cursor();
  h->ng_probs->init_cursor();
  while(ng_lst->next()){
    h->ng_probs->next();
    h->ng_probs->write(ng_lst->read()->get_prob(h,nextword));
  }
}


DOUBLE
MIXTURE_NGRAM::text_score(QUEUE<INT> *wq, GENFILE *outfp, LEXICON *lex, 
		  DOUBLE &noUNK_accum, INT &noov) {
  DOUBLE score = 0.0;
  LIST<INT> *cxt = new LIST<INT>();
  QUEUE<INT> *newq = new QUEUE<INT>();
  INT wd;
  DOUBLE oldp = noUNK_accum;
  HYP *hy;


  if(Debug_flag.text_decode)
    fprintf(stderr, "In text_score()\n");

  // bbo default cxt is </s>
  // but for noway it is <s>
  for(INT j = 0; j < order-1; j++)
    cxt->push(Sentence_start_index);
    //    cxt->push(Sentence_end_index);

  HYP *rooth = new HYP(Sentence_start_index, 0, 0, 0, NULL, 
		       oov2unk(Sentence_start_index));
  hy = rooth;
 
  while(wq->size > 0) {
    wd = wq->pop();
    if(wd != Sentence_start_index && wd != Paragraph_index && wd != Article_index) {
      //      INT p = get_prob(cxt, wd);
      INT p = get_prob(hy, wd);
      //      score += p;
      score += logdecodelog(p);
      if(Verbose) {
	printf("log_10 P(%s | ", lex->vocabulary->get(wd));
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
    hy = new HYP(wd, 0, 0, 0, hy, oov2unk(wd));
  }
  cxt->clear();
  outfp->fprintf(" (%.2f, %.2f including <UNK>)\n", noUNK_accum - oldp, 
	  score);
	  //	  logdecodelog(score)); 
  delete newq;
  delete cxt;
  return score;
}
