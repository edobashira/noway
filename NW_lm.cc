// -*- Mode: C++;  -*-
// File: lm.cc
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
#include <math.h>
#include <string.h>

// Parameterised classes
#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_debug.h"
#include "NW_param.h"

#include "NW_hypothesis.h"
#include "NW_lm.h"
#include "NW_lexicon.h"


BOOL 
GRAMMAR::hyp_exact_similar(LIST<INT> *h1h, LIST<INT> *h2h) {
  INT w1, w2;

  h1h->init_cursor();
  h2h->init_cursor();
  while(h1h->next() && h2h->next()) {
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
    if(w1 == -1)
      break;
  }
  if(h1h->end_cursor() && h2h->end_cursor())
    return(TRUE);
  
  return(FALSE);
}


// Doesn't work with online_output
BOOL
GRAMMAR::hyp_exact_similar(HYP *h1, HYP *h2){
  while(h1 != NULL && h2 != NULL) {
    if(h1->get_word() == Silence_index)
      h1 = h1->get_previous();
    if(h2->get_word() == Silence_index)
      h2 = h2->get_previous();
    if((h1 == NULL || h2 == NULL) ||
       h1->get_word() != h2->get_word())
      break;
    h1 = h1->get_previous();
    h2 = h2->get_previous();
  }
  return(h1 == NULL && h2 == NULL);
}


// Doesn't work with online_output
BOOL 
GRAMMAR::hyp_exact_similar(LIST<INT> *h1h, LIST<INT> *h2h, INT w) {
  INT w1, w2;
  
  h1h->init_cursor();
  h2h->init_cursor();
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

  while(h1h->next() && h2h->next()) {
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
    if(w1 == -1)
      break;
  }
  if(h1h->end_cursor() && h2h->end_cursor())
    return(TRUE);

  return(FALSE);
}

BOOL
GRAMMAR::hyp_exact_similar(HYP *h1, HYP *h2, INT w){
  HYP* tmph = new HYP(w, 0, 0, 0, h2);
  BOOL res = hyp_exact_similar(h1, tmph);
  delete tmph;
  return(res);
}



DOUBLE
GRAMMAR::text_decode(GENFILE *fp, GENFILE *outfp, LEXICON *lex, 
		     DOUBLE &noUNKp, BOOL ignore_newline) {
  CHAR buf[BUFSIZ];
  CHAR *cptr, *wptr;;
  QUEUE<INT> *wq = new QUEUE<INT>();
  INT word_id = -1;
  DOUBLE perplexity;
  DOUBLE noUNK_perplexity;
  DOUBLE log_pp = 0.0;
  DOUBLE noUNK_log_pp = 0.0;
  //  INT prob_accum = 0;
  DOUBLE prob_accum = 0.0;
  INT nwords = 0;
  INT nsents = 0;
  INT ntags = 0;
  //  INT noUNK_prob_accum = 0;
  DOUBLE noUNK_prob_accum = 0.0;
  INT noov = 0;
  INT nw_adjust = 0;
  
  if(use_unk)
    outfp->fprintf("unk_factor = %g\n", 1.0/logdecode(unk_factor));

  // ignore newline
  if(ignore_newline) {
    while(!fp->eof()) {
      wq->fast_clear();
      do {
	if(fscanf(fp->get_fp(), "%s", buf) != 1)
	  break;
	word_id = lex->vocab_hash->get(buf);
	if(word_id == -1) {
	  if(use_unk)
	    word_id = Unknown_word_index;
	  else
	    panic("text_decode(), closed vocab, unknown word %s", buf);
	}
	// Do we want this??
	if(wq->size == 0 && (word_id != Sentence_start_index || word_id != Paragraph_index || word_id != Article_index))
	  wq->push(Sentence_start_index);
	wq->push(word_id);	
      } while(word_id != Sentence_end_index);

      if(wq->size == 0) break;

      if(word_id != Sentence_end_index){
	wq->push(Sentence_end_index);
	nw_adjust++;
      }

      nsents++;
      nwords += wq->size;
      //      noUNK_prob_accum = 0;
      noUNK_prob_accum = 0.0;
      prob_accum = text_score(wq, outfp, lex, noUNK_prob_accum, noov);


      //      log_pp += logdecodelog(prob_accum);
      //      noUNK_log_pp += logdecodelog(noUNK_prob_accum);
      log_pp += prob_accum;
      noUNK_log_pp += noUNK_prob_accum;
    }
  } else {
    // treat newline as </s>
    while(fp->fgets(buf, BUFSIZ) != NULL) {
      wq->fast_clear();
      wq->push(Sentence_start_index);
      // miss out comments
      if(buf[0] == '#')
	continue;

      nsents++;
      //      fprintf(outfp, "%s", buf);
      cptr = buf;
      while((wptr = strtok(cptr, " \t\n")) != NULL) {
	cptr = NULL;
	word_id = lex->vocab_hash->get(wptr);
	if(word_id == -1) {
	  if(use_unk)
	    word_id = Unknown_word_index;
	  else
	    fprintf(stderr, "text_decode(), closed vocab, unknown word %s", wptr);
	}
	wq->push(word_id);
      }
      if(word_id != Sentence_end_index){
	wq->push(Sentence_end_index);
	nw_adjust++;
      }

      nwords += wq->size;
      //      noUNK_prob_accum = 0;
      noUNK_prob_accum = 0.0;
      prob_accum = text_score(wq, outfp, lex, noUNK_prob_accum, noov);
      //      log_pp += logdecodelog(prob_accum);
      //      noUNK_log_pp += logdecodelog(noUNK_prob_accum);
      log_pp += prob_accum;
      noUNK_log_pp += noUNK_prob_accum;
    }
  }

  // subtract off <s> perplexity calcn.
  nwords -= (nsents+nw_adjust);

  perplexity = exp10(-log_pp/(DOUBLE)nwords);
  INT noUNK_nwords = nwords - noov;
  if(noUNK_nwords > 0)
    noUNK_perplexity = exp10(-noUNK_log_pp/(DOUBLE)noUNK_nwords);
  else
    noUNK_perplexity = 0.0;

  outfp->fprintf("sum_logprob = %.2f (%.2f including <UNK>)\n", noUNK_log_pp, log_pp);
  outfp->fprintf("%d probs evaluated (%d including <UNK>)\n", noUNK_nwords, nwords);
  outfp->fprintf("\nTest Set Perplexity = %.2f (%.2f including <UNK>)\n", 
	  noUNK_perplexity, perplexity);

  // subtract off </s> for OOV calcn.
  //nwords -= nsents;
  //  noUNK_nwords -= nsents;

  INT nreal_words = nwords-ntags;
  outfp->fprintf("%d words total:  %d OOV (%.2f%%)\n", nwords, noov, 
		 (FLOAT)noov/(FLOAT)nreal_words);
  outfp->fprintf("                 %d IV words\n", nreal_words-noov);  
  noUNKp = noUNK_perplexity;
  return perplexity;
}

// outfp->fprintf("%d out of %d words (%.2f%%) OOV\n", nwords - noUNK_nwords,
//	  nwords, 100.0 * (DOUBLE)(nwords - noUNK_nwords)/((DOUBLE)nwords));
