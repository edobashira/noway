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
#include "NW_tagged_ngram.h"
#include "NW_lexicon.h"

// \LOCATION\
// ngram 1 = 1236
//  
// \1-grams:
// -3.8130 <UNK>
// [...]
// -3.5116 YUCCA
// 
// \end of LOCATION\
// 
// \NAME\
// [...]



INT
TAG_UNIGRAM::build(GENFILE *fp, LEXICON *lex) {
  char buf[BUFSIZ];
  char *tagname;
  INT indx = -1;
  INT nw = lex->vocabulary->size;
  INT nread;
  LIST<INT>* probs;
  FLOAT prob;


  if(Verbose)
    printf("Building TAG_UNIGRAM\n");
  nread = fscanf(fp->get_fp(), " %s ", buf);
  while(nread == 1 && !fp->eof()){
    if(buf[0] != '\\' && buf[strlen(buf)-1] != '\\')
      panic("TAG_UNIGRAM::build() - Failed reading %s\n", buf);
    buf[strlen(buf)-1] = '>';
    buf[0] = '<';
    tagname = strdup(buf);
    INT order, nprobs;
    probs = new LIST<INT>(nw);
    for(INT j = 0; j < nw; j++)
      probs->push(PROB_NULL);
    tag_probs->push(probs);
    tag_strs->push(tagname);
    indx =lex->vocab_hash->get(tagname);
    if(indx < 0) 
      panic("Tag %s not in ngram model!", tagname);
    tag_ids->push(indx);

    nread = fscanf(fp->get_fp(), " ngram %d = %d", &order, &nprobs);
    if(nread != 2)
      panic("TAG_UNIGRAM::build() - Failed to read ngram field for tag %s", tagname);
    nread = fscanf(fp->get_fp(), " ngram %d = %d", &order, &nprobs);
    if(nread != 0)
      panic("TAG_UNIGRAM::build() - Can only deal with uningram tag probs");

    nread = fscanf(fp->get_fp(), " \\%s", buf);
    if(nread != 1 || strcmp(buf, "1-grams:"))
      panic("TAG_UNIGRAM::build() - Failed, expecting \\1-grams: line");
    if(Verbose)
      printf("Processing %d %s probs\n", nprobs, tagname);
    for(INT i = 0; i < nprobs; i++){
      nread = fscanf(fp->get_fp(), " %f %s", &prob, buf);
      if(nread != 2)
	panic("TAG_UNIGRAM::build() - Failed to read %dth prob for tag %s\n",
	      i, tagname);
      indx = lex->vocab_hash->get(buf);
      if(indx == -1) {
	indx = lex->vocabulary->size | LM_OOV_MASK;
	lex->vocabulary->push(strdup(buf));
	lex->vocab_hash->insert(strdup(buf), indx);
      } 
      probs->put_and_grow(indx&WORD_MASK, MAX(prob, -7.999) * log_scale + 0.5, 
			    PROB_NULL);
    }
    nread = fscanf(fp->get_fp(), " \\end of %s", buf+1);
    if(nread != 1)
       panic("TAG_UNIGRAM::build() - Failed to read end of %s", tagname);
    buf[0] = '<';
    buf[strlen(buf)-1] = '>';
    if(strcmp(buf, tagname))
      panic("TAG_UNIGRAM::build() - Mismatch \\%s\\ \\end of %s\\", tagname, buf);
    delete tagname;
    nread = fscanf(fp->get_fp(), " %s ", buf);
  }
  if(Verbose)
    printf("Built TAG_UNIGRAM, %d tags\n", tag_strs->size);
  
  return(lex->vocabulary->size);
}

VOID
TAG_UNIGRAM::info(GENFILE* f) {
}


INT 
TAGGED_NGRAM::get_prob(HYP* h, INT nextword, INT hypid){
  INT res = 0;
  INT p;
  INT t;
  INT indx = nextword & WORD_MASK;

   
  if(nextword & LM_OOV_MASK){
    res = logcodelog(-99.9);
    tag->prob_init_cursor();
    tag->id_init_cursor();
    while(tag->prob_next()){
      tag->id_next();
      if(indx < tag->prob_read()->size) {
	p = tag->prob_read()->get(indx);
	if(p != PROB_NULL){
	  t = tag->id_read();
	  res = logadd(res, p + ng->get_prob(h, t, hypid));
	}
      }
    }
  } else {
    // nextword is not OOV
    res = ng->get_prob(h, nextword, hypid);    
  }
	
  h->lm_prob = res;
  return res;  
}

INT 
TAGGED_NGRAM::get_prob(LIST<INT>* wlist, INT nextword, INT hypid){
  INT res = 0;
  INT p;
  INT t;
  INT indx = nextword & WORD_MASK;

   
  if(nextword & LM_OOV_MASK){
    res = logcodelog(-99.9);
    tag->prob_init_cursor();
    tag->id_init_cursor();
    while(tag->prob_next()){
      tag->id_next();
      if(indx < tag->prob_read()->size) {
	p = tag->prob_read()->get(indx);
	if(p != PROB_NULL){
	  t = tag->id_read();
	  res = logadd(res, p + ng->get_prob(wlist, t, hypid));
	}
      }
    }
  } else {
    // nextword is not OOV
    res = ng->get_prob(wlist, nextword, hypid);    
  }
	
  return res;  
}

DOUBLE 
TAGGED_NGRAM::text_score(QUEUE<INT> *wq, GENFILE *outfp, LEXICON *lex, 
			 DOUBLE &noUNK_accum, INT &noUNK_nwords) {
  DOUBLE score = 0.0;
  LIST<INT> *cxt = new LIST<INT>();
  INT wd;
  DOUBLE oldp = noUNK_accum;

  if(Debug_flag.text_decode)
    fprintf(stderr, "In text_score()\n");

  // bbo default cxt is </s>
  // but for noway it is <s>
  for(INT j = 0; j < ng->get_order()-1; j++)
    cxt->push(Sentence_start_index);
    //    cxt->push(Sentence_end_index);

  while(wq->size > 0) {
    wd = wq->pop();
    if(wd != Sentence_start_index && wd != Paragraph_index && wd != Article_index) {
      INT p, tp, t; 
      //      score += p;
      INT indx = wd&WORD_MASK;
      if(Verbose) {
	printf("P(%s | ", lex->vocabulary->get(indx));
	cxt->init_cursor();
	for(INT count = 1; count < ng->get_order() && cxt->next(); count++)
	  printf("%s ", lex->vocabulary->get(cxt->read()));
	printf(") = ");
	if(wd & LM_OOV_MASK){
	  p = logcodelog(-99.9);
	  tag->prob_init_cursor();
	  tag->id_init_cursor();
	  tag->strs_init_cursor();
	  while(tag->prob_next()){
	    tag->id_next();
	    tag->strs_next();
	    if(indx < tag->prob_read()->size) {
	      tp = tag->prob_read()->get(indx);
	      if(tp != PROB_NULL){
		t = tag->id_read();
		p = logadd(p, tp + ng->get_prob(cxt, t));
		printf("P(%s | %s) = %f  ", lex->vocabulary->get(indx), 
		       tag->strs_read(), 
		       logdecodelog(tp));
		printf("P(%s | ", tag->strs_read());
		cxt->init_cursor();
		for(INT count = 1; count < ng->get_order() && cxt->next(); count++)
		  printf("%s ", lex->vocabulary->get(cxt->read()));
		printf(") = %f", logdecodelog(tp + ng->get_prob(cxt, t)));
	      }
	    }
	  }
	  printf(" = ");
	} else {
	  p = ng->get_prob(cxt, wd);
	  printf("P(%s | ", lex->vocabulary->get(indx));
	  cxt->init_cursor();
	  for(INT count = 1; count < ng->get_order() && cxt->next(); count++)
	    printf("%s ", lex->vocabulary->get(cxt->read()));
	  printf(") = ");
	}
	printf("%f\n", logdecodelog(p));
      } else {
	p = ng->get_prob(cxt, wd);
      }
      score += logdecodelog(p);
      if(wd != Unknown_word_index) {
	noUNK_accum += logdecodelog(p);
	noUNK_nwords++;
      }
    }
    cxt->push(oov2unk(wd));
  }
  cxt->clear();
  outfp->fprintf(" (%.2f, %.2f including <UNK>)\n", noUNK_accum - oldp, 
	  score);
  delete cxt;
  return score;
}










