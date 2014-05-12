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



// Parameterised classes
#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_debug.h"
#include "NW_param.h"

#include "NW_hypothesis.h"

void* 
HYP::operator new(size_t size) {
  if(size != sizeof(HYP))
    return ::new char[size];

  HYP* res = free_list;

  if(res == NULL){
    // add more memory to the free list

    // use ::new to avoid calling the constructor
    HYP* new_block = 
      (HYP*) ::new char[block_size * sizeof(HYP)];

    if(new_block == NULL) return NULL;
    
    for(INT i = 1; i < block_size; i++){
      new_block[i-1].previous = &new_block[i];
      // for debugging
      //      new_block[i].n_next = -100;
    }
    new_block[block_size-1].previous = NULL;

    for(INT j = 0; j < block_size; j++){
      new_block[j].lm_context = new INT[lm_context_length];
      if(State_decode) {
	new_block[j].nodeseq = new QUEUE<INT>();
	new_block[j].nodetimeseq = new QUEUE<INT>();
	new_block[j].nodeprobseq = new QUEUE<INT>();
      } else {
	new_block[j].nodeseq = NULL;
	new_block[j].nodetimeseq = NULL;
	new_block[j].nodeprobseq = NULL;
      }
      if(mixture_lm_size > 0){
	//	new_block[j].next_lm_wts = new LIST<INT>(mixture_lm_size);
	//	new_block[j].current_lm_wts = new LIST<INT>(mixture_lm_size);
	new_block[j].ng_wts = new LIST<INT>(mixture_lm_size);
      } else {
	//	new_block[j].next_lm_wts = NULL;
	//	new_block[j].current_lm_wts = NULL;
	new_block[j].ng_wts = NULL;
      }
      if(gamma == NULL) {
	gamma = new LIST<INT>(mixture_lm_size);
	gamma->set(0,mixture_lm_size);
      }
      if(ng_probs == NULL) {
	ng_probs = new LIST<INT>(mixture_lm_size);
	ng_probs->set(0,mixture_lm_size);
      }
    }

    res = new_block;

  }
  
  free_list = res->previous;
  return res;
}


void
HYP::operator delete(void* dead, size_t size){
  if(size != sizeof(HYP)){
    ::delete [] ((char*) dead);
    return;
  } 
  ((HYP*)dead)->previous = free_list;
  free_list = (HYP*)dead;
}


// EM blind re-estimation of c^t_j:
// c^t_j = 1/t \sum_{u=1}^t \gamma_j^u
// \gamma_j^t = \frac{c_j^{t-1} tg(w | w_h; M_j)}
//                   {sum_k c_k^{t-1} tg(w | w_h; M_k)}
// This assumes that HYP::ng_probs->get(j) is set to tg(w | w_h; M_j)
//  and that HYP::lm_prob is set to the wtd sum
VOID
HYP::update_ng_wts() {
  INT logcount = 0;
  INT logratio = -logadd_lookup_length;
  LIST<INT>* oldwts = NULL;

  assert(previous != NULL);
  oldwts = previous->get_ng_wts();
 
  assert(mixture_lm_size > 0);
  assert(ng_probs->size == mixture_lm_size);
  assert(oldwts->size == mixture_lm_size);
  assert(gamma->size == mixture_lm_size);

  if(count > 2){
    logcount = logcode((FLOAT)(count-1));
    logratio = logcode((FLOAT)count-2) - logcount;
  }
    
  gamma->init_cursor();
  ng_probs->init_cursor();
  while(gamma->next()){
    ng_probs->next();
    gamma->write(ng_probs->read() - lm_prob);
  }

  if(Debug_flag.mixlm) {
    printf("\n\ncount = %d\n", count);
    gamma->init_cursor();
    ng_probs->init_cursor();
    while(gamma->next()){
      ng_probs->next();
      printf("  gamma = %f,  loggamma = %f\n",
	     logdecode(gamma->read()), logdecodelog(gamma->read()));
    }
  }
  ng_wts->init_cursor();
  gamma->init_cursor();
  oldwts->init_cursor();
  while(ng_wts->next()){
    gamma->next();
    oldwts->next();
    ng_wts->write(logadd(oldwts->read()+logratio, gamma->read()-logcount));
  }

  if(Debug_flag.mixlm && count==3){
    ng_wts->set(logcode(0.1), mixture_lm_size);
    printf("ng_wts->set(%d/%f, %d)\n", logcode(0.1), 0.1, mixture_lm_size);
  }
}

// This assumes that ng_probs and ng_wts have been set
// sets lm_prob
INT
HYP::wt_sum_probs() {
  INT x;
  assert(ng_probs->size == mixture_lm_size);
  assert(ng_wts->size == mixture_lm_size);
  assert(mixture_lm_size > 0);
  INT res = ng_probs->get(0) + ng_wts->get(0) - logadd_lookup_length;
  ng_probs->init_cursor();
  ng_wts->init_cursor();
  while(ng_probs->next()){
    ng_wts->next();
    if(Debug_flag.mixlm)
      printf("logFactor = %f Factor %f   logP = %f P = %f\n", 
	     logdecodelog(ng_wts->read()), logdecode(ng_wts->read()),
	     logdecodelog(ng_probs->read()), logdecode(ng_probs->read()));
    x = ng_probs->read() + ng_wts->read();
    res = logadd(res, x);
    ng_probs->write(x);
  }
  if(Debug_flag.mixlm)
      printf("Total logP = %f P = %f\n", logdecodelog(res), logdecode(res));
  lm_prob = res;
  return(res);  
}



VOID
HYP::display(LIST<STRING> *vocab, GENFILE *fp, INT sc,
		     BOOL op_logprob, STRING first_char){
  INT w;
  LIST<INT>* wlist = new LIST<INT>();
  LIST<INT>* vlist = new LIST<INT>();
  HYP* p;

  fp->fprintf("%s", first_char);
  if(sc != -1)
    fp->fprintf("%d ", sc);

  // build the word list - first word will end up on top
  p = this;
  while(p != NULL) {
    w = p->get_word();
    if(w != Silence_index && w != Sentence_start_index && w != Sentence_end_index){
      wlist->push(w&WORD_MASK);
      if(output_pron_version)
	vlist->push(p->get_pron_version());
    }
    p = p->get_previous();
  }
  
  while(wlist->size > 0){
    if(output_pron_version)
      fp->fprintf("%s(%d) ", vocab->get(wlist->pop()), vlist->pop());
    else
      fp->fprintf("%s ", vocab->get(wlist->pop()));
  }
  delete wlist;
  delete vlist;

  if(op_logprob)
    fp->fprintf(" (%f)", logdecodelog(logprob));
  fp->fprintf("\n");
  fp->fflush();
}

VOID
HYP::log_display(LIST<STRING> *vocab, GENFILE *fp){
  INT w;
  LIST<INT>* wlist = new LIST<INT>();
  LIST<INT>* vlist = new LIST<INT>();
  LIST<INT>* tlist = new LIST<INT>();
  LIST<INT>* plist = new LIST<INT>();
  HYP* p;

  // build the word list - first word will end up on top
  p = this;
  while(p != NULL) {
    w = p->get_word();
    wlist->push(w&WORD_MASK);
    vlist->push(p->get_pron_version());
    plist->push(p->get_logprob());
    tlist->push(p->get_time());
    p = p->get_previous();
  }
  
  while(wlist->size > 0){
    w = wlist->pop();
    if(vocab)
      fp->fprintf("%s/%d(%d) [%.2f] <%d> ", vocab->get(w), w, 
		  vlist->pop(),
		  logdecodelog(plist->pop()), tlist->pop());
    else
      fp->fprintf("%d(%d) [%.2f] <%d> ", w, vlist->pop(), 
		  logdecodelog(plist->pop()), tlist->pop());
  }
  delete wlist;
  delete plist;
  delete tlist;
  delete vlist;
  fp->fprintf("\n");
  fp->fprintf("  time = %d, logprob = %d (%f), deltalogprob = %d\n", get_time(), 
	      get_logprob(), logdecodelog(get_logprob()), get_delta());
  fp->fflush();
}

VOID
HYP::state_display(LIST<STRING> *phnms, INT offset, FLOAT fshift, GENFILE *fp){
  INT w, start, end;
  LIST<INT>* statel = new LIST<INT>();
  LIST<INT>* timel = new LIST<INT>();
  HYP* p;

  // build the state list - first state will end up on top
  p = this;
  while(p != NULL) {
    p->nodeseq->init_cursor();
    while(p->nodeseq->next())
      statel->push(p->nodeseq->read());
    p->nodetimeseq->init_cursor();
    while(p->nodetimeseq->next())
      timel->push(p->nodetimeseq->read());
    p = p->get_previous();
  }

  assert(timel->size == statel->size);
  
  start = 0;
  end = 0;
  while(statel->size > 0) {
    w = statel->pop();
    end = timel->pop();
    if(start == end)
      continue;
    fp->fprintf("%.3f ", (start+offset) * fshift);
    fp->fprintf("%.3f ", (end-start) * fshift);
    fp->fprintf("%s (%d)\n", phnms->get(w), w);
    start = end;
  }
  delete statel;
  delete timel;
  fp->fprintf("\n");
  fp->fflush();  
}

VOID
HYP::ctm_display(LIST<STRING> *vocab, INT offset, FLOAT fshift,
			 GENFILE *fp, INT sc, CHAR channel){
  INT w, pr, start, end, v;
  LIST<INT>* wlist = new LIST<INT>();
  LIST<INT>* tlist = new LIST<INT>();
  LIST<INT>* plist = new LIST<INT>();
  LIST<INT>* vlist = new LIST<INT>();
  HYP* p;
  
  if(Pause_decode == FALSE)
    return;

  // build the word list - first word will end up on top
  p = this;
  while(p != NULL) {
    w = p->get_word();
    if(w == Sentence_start_index || w == Sentence_end_index)
      w = Silence_index;
    wlist->push(w&WORD_MASK);
    plist->push(p->get_word_logprob());
    tlist->push(p->get_time());
    vlist->push(p->get_pron_version());
    p = p->get_previous();
  }

  start = 0;
  while(wlist->size > 0) {
    end = tlist->pop();
    pr = plist->pop();
    v = vlist->pop();
    w = wlist->pop();
    if(w == (Silence_index&WORD_MASK) && !output_silence)
      continue;
    if(start == end)
      continue;
    if(channel == '#')
      fp->fprintf("%c ", channel);
    else
      fp->fprintf("%d %c ", sc, channel);
    fp->fprintf("%.3f ", (start+offset) * fshift);
    fp->fprintf("%.3f ", (end-start) * fshift);
    if(output_pron_version)
      fp->fprintf("%s (%d) ", vocab->get(w), v);
    else
      fp->fprintf("%s ", vocab->get(w));
    fp->fprintf("%.4f", logdecodelog(pr));
    fp->fprintf("\n");   
    start = end;
  }
  fp->fflush();
  delete wlist;
  delete plist;
  delete tlist;
  delete vlist;
}

VOID
HYP::srt_display(LIST<STRING> *vocab, INT offset, FLOAT fshift, GENFILE *fp){
  INT w, pr, start, end, v;
  LIST<INT>* wlist = new LIST<INT>();
  LIST<INT>* vlist = new LIST<INT>();
  LIST<INT>* tlist = new LIST<INT>();
  LIST<INT>* plist = new LIST<INT>();
  HYP* p;
  
  if(Pause_decode == FALSE)
    return;

  // build the word list - first word will end up on top
  p = this;
  while(p != NULL) {
    w = p->get_word();
    if(w == Sentence_start_index || w == Sentence_end_index)
      w = Silence_index;
    wlist->push(w&WORD_MASK);
    plist->push(p->get_word_logprob());
    tlist->push(p->get_time());
    vlist->push(p->get_pron_version());
    p = p->get_previous();
  }

  start = 0;
  while(wlist->size > 0) {
    w = wlist->pop();
    end = tlist->pop();
    pr = plist->pop();
    v = vlist->pop();
    if(w == (Silence_index&WORD_MASK) && !output_silence)
      continue;
    if(end-start == 0)
      continue;
    fp->fprintf("<Word ");
    fp->fprintf("S_time=%.3f ", (start+offset) * fshift);
    fp->fprintf("E_time=%.3f ", (end+offset) * fshift);
    fp->fprintf("Prob=%.4g> ", logdecodelog(pr));
    if(output_pron_version)
      fp->fprintf("Version=%d> ", v);
    fp->fprintf("%s </Word>", vocab->get(w));
    fp->fprintf("\n");
    start = end;
  }
  fp->fflush();  
  delete wlist;
  delete plist;
  delete tlist;
  delete vlist;
}


VOID
HYP::info(LIST<STRING> *vocab, LIST<STRING> *phns, 
		  GENFILE *fp, STRING hdr){
  fp->fprintf("%s\n", hdr);
  fp->fprintf("Hypothesis element info (%d total hyps)\n", get_total_hyps());
  fp->fprintf("    address = %p,   ", this);
  fp->fprintf("   previous = %p\n", previous);
  fp->fprintf("   time = %d, ", t);
  if(vocab) 
    fp->fprintf("   word = %s (%d), version %d\n", vocab->get(wd), wd, pron_version);
  else
    fp->fprintf("   word = %d, version %d\n", wd, pron_version);
  fp->fprintf("   logprob = %g (%d), ", logdecodelog(logprob), logprob);
  fp->fprintf("   wd_logprob = %g (%d), ", logdecodelog(wd_logprob), wd_logprob);
  fp->fprintf("   pron_prob = %g (%d), ", logdecodelog(pron_prob), pron_prob);
  fp->fprintf("   deltalogprob = %d\n", deltalogprob);
  fp->fprintf("   n_next = %d,  ", get_n_next());
  fp->fprintf("LM context: ");
  if(vocab) {
    for(INT i = 0; i < lm_context_length; i++)
      fp->fprintf(" %s/%d ", vocab->get(get_lm_context()[i]), get_lm_context()[i]);
  } else {
    for(INT i = 0; i < lm_context_length; i++)
      fp->fprintf(" %d ", get_lm_context()[i]);
  }
  HYP* hp = get_previous();
  INT np = 0;
  while(hp != NULL){
    np++;
    hp = hp->get_previous();
  }
  fp->fprintf("%d previous in chain\n", np);  
  //  if(vocab){
  //    fp->fprintf("Backtracing:\n");
  //    display(vocab, fp);
  //  }
  if(phns)
    if(State_decode)
      state_display(phns, 0, 1., fp);
  fp->fprintf("\n");
}
