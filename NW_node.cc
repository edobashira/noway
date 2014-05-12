// -*- Mode: C++;  -*-
// File: node.cc
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Forward, prune, create_new_hyps functions of node
//*             in lexical tree
//*
//* CLASSES:
//* 
//* REQUIRES:
//*
//* HISTORY:
//*  Jul 20 17:25 1994 (sjr): added fast lookahead (checking posterior
//*                           floor for min duration before activating
//*                           a node)
//*  May 24 15:46 1994 (sjr): Test okay in decoder
//* Created: Thu Apr 14 10:59:44 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "all.h"

NODE*
NODE::match_succ(HMM *p) {
  NODE* res;
  for(INT i = succ->size-1; i >= 0; i--){
    res = succ->get(i);
    if(res->is_equal(p)) 
      return(res);
  }
  return(NULL);
}

NODE*
NODE::create_successor(HMM *p, INT nh, BOOL hp) {
  NODE* res = new NODE(p, this, nh, hp);
  succ->push(res);
  return(res);
}


// Given a list of hypotheses and a time t, this function makes one 
//  step forward in time, computing new node occupation probs,
//  activating new nodes if necessary and registering any word
//  endings.  Finally a call is made to the pruner.
// When a node is deactivated, the log_prob array is set to PROB_NULL = 0x7f7f7f7f
//  log_prob[0] is the entry probability at time t for this node
//  log_prob[1] is the exit probability at time t for this node
// 
// Returns TRUE if node is to be entered on the active queue at the next time
// called from TREE_LEXICON::extend() (lex_search.cc)
BOOL 
NODE::forward(INT t, LIST<HYP*>* hyps, GRAMMAR* lm,
	      QUEUE<QUEUE<NODE*>*>* act, INT& maxlevel) {
  // step_forward() returns FALSE when out of data (so don't consider 
  //  successor nodes)
  if(!step_forward(t))
    return(FALSE);

  // If exit state has a probability, then get the list of successors
  //  For each successor node:
  //   if(not active) activate, and add to active queue act
  //   else(active) set logprob[ENTRY_STATE] and add to nodeq, which 
  //    will be appended to the active queue for the next time step
  //  NB breadth-first queue ensures parent always precedes children
  INT lp = logprob[EXIT_STATE];
  if(lp != PROB_NULL)  {
    NODE* nd;
    INT op;
    INT lmapprox = 0;
    for(INT i = 0; i < succ->size; i++){
      nd = succ->get(i);
      op = nd->output_prob(t);
      // look ahead for floor or end of sentence
      // note this assumes strict left-right models with no skips...
      if(is_floor(op) || t+nd->get_nstates()-FIRST_STATE >= sentence_length)
	continue;

      if(nd->retrieved_lmprob == TRUE) {
	lmapprox = nd->maxlmub;
      } else if (retrieved_lmprob == TRUE) {
	lmapprox = maxlmub;
      } else {
	lmapprox = 0;
	if(smear_unigram == TRUE)
	  lmapprox += nd->maxuniprob;
	if(smear_context == TRUE)
	  lmapprox += maxlmub;  // maxlmub is same for all nodes
      }

      if(delta(op + lp + lmapprox, t) < state_beam) {
	if(nd->activated){
	  // Each node only has one parent, and cannot reactivate 
	  //  its own entry state
	  nd->set_init_logprob(lp, t);
	} else {
	  BOOL addnd = TRUE;
	  INT j, k;
	  for(k = FIRST_STATE+1, j = t+1; k < nd->get_nstates(); k++, j++) {
	    if(is_floor(nd->output_prob(j))) {
	      addnd = FALSE;
	      break;
	    }
	  }
   	  if(addnd){
	    nd->activate(t);
 	    //	    if(act->size == 0 || nd->level >= act->last()->level) {
	    //	      act->enqueue(nd);
	    //	    } else {
	    //	      nreactivated++;
	    //	      excep->push(nd);	      
	    //	    }
	    act->get(nd->level)->push(nd);
	    maxlevel = MAX(maxlevel, nd->level);
	  }
	}
      }
    }
  }

  // check for node completion and pruning
  return(!prune(t, hyps, lm));
}

//  This prunes nodes, relative to hypotheses, that are out of the likelihood 
//   beam, using the grammar.  If a node is pruned realtive to all 
//   hypotheses then that node is pruned and TRUE is returned.
BOOL
NODE::prune(INT t, LIST<HYP*> *hyps, GRAMMAR *lm) {
  INT j, p;
  HYP *hy;
  INT lmp;
  INT lmapprox = 0;
  assert(hyps->size > 0);

  // Compute the upper bound LM probability for each hypothesis at this node
  //  ie for each hypothesis the LM probability of the most probable
  //  pathword.  If too many pathwords, ignore the LM in this node.
  if(pathword_threshold > 0 && 
     pathword->size <= pathword_threshold && 
     retrieved_lmprob == FALSE){
    assert(maxlmprob != NULL);
    maxlmub = INT_MIN;
    maxlmprob->set(INT_MIN, hyps->size);
    hyps->init_cursor();
    if(hyp_active) {
      assert(hyp_active->size == hyps->size);
      hyp_active->init_cursor();
    }
    maxlmprob->init_cursor();
    INT hcount = -1;
    while(hyps->next()){
      hcount++;
      if(hyp_active) hyp_active->next();
      maxlmprob->next();
      // this relies on left-right evaluation
      if(hyp_active == NULL || hyp_active->read()) {
	hy = hyps->read();
	for(j = 0; j < pathword->size; j++){
	  p = lm->get_prob(hy, pathword->get(j), hcount) - hy->get_delta();
	  if(p > maxlmprob->read()) {
	    maxlmprob->write(p);
	    maxlmub = MAX(p, maxlmub);
	  }	 
	}
      }
    }
    retrieved_lmprob = TRUE;
    lmapprox = maxlmub;
  } else {
    lmapprox = 0;
    if(smear_unigram == TRUE)
      lmapprox += maxuniprob;
    if(smear_context == TRUE)
      lmapprox += maxlmub;
  }

  BOOL pru = TRUE;
  INT d, mind = INT_MAX;

  // First try using lmapprox/maxlmub to see if the node is pruned
  maxprob = INT_MIN;
  for(j = nstates-1; j >= FIRST_STATE; j--){
    if(logprob[j] == PROB_NULL) 
      continue;
    p = logprob[j] + lmapprox;
    d = delta(p, t);
    if(d <= state_beam) {
      pru = FALSE;
      if(!strict_lub && maxprob < p) {
	maxprob = p;
	mind = d;
      } else if(logprob[j] > Bestprob_this_ext) {
	Bestprob_this_ext = logprob[j];
	Best_t = t;
      }
    } else {
      logprob[j] = PROB_NULL;
    }
  }
  // greedy lub update
  if(!strict_lub && mind < 0) 
    update_lub(maxprob, t);

  if(Debug_flag.prune)
    printf("%d(t=%d): maxp = %d, lub = %d, prune = %d\n", 
	   id, t, maxprob, lub->get(t), pru);


  if(pru == FALSE && retrieved_lmprob == TRUE && hyp_active != NULL) {
    // Now prune individual hypotheses 
    INT lowest_unpruned = maxlmub;
    assert(hyp_active->size == maxlmprob->size);
    for(INT i = 0; i < maxlmprob->size; i++){
      lmp = maxlmprob->get(i);
      if(hyp_active->get(i) && lmp < lowest_unpruned) {
	hyp_active->put(i, FALSE);
	if(delta(maxprob+lmp-maxlmub, t) < state_beam) {
	  lowest_unpruned = lmp;
	  hyp_active->put(i, TRUE);
	  break;
	}
      }
    }
  }

  return(pru);
}

INT
NODE::create_new_hyps(INT t, LIST<HYP*> *hyps, GRAMMAR *lm,
		      QUEUE<STACK*> *stks) {
  HYP* h = NULL;
  HYP* newh = NULL;
  STACK *s;
  INT p, lmp, wp, ap, w, offp;
  INT v; // DCA 26/FEB/96: pronunciation version variable
  INT res = 0;
  INT nadded = 0;

  // Should not be called if exit state is not active or node has no endwords
  assert(!(logprob[EXIT_STATE] == PROB_NULL || endword->size == 0 || t == 0));

  INT index = t - (start_time + 1);
  if(index < 0) {
    if(t == start_time)
      return(0);
    else
      panic("t = %d, start_time = %d: Rogue hypothesis in node %s(%d)\n", 
	    t, start_time, name, id);
  }

  // create enough new stacks to reach the ref time of the new hypotheses
  while(index >= stks->size)
    stks->enqueue(new STACK());

  // get the stack with the relevant ref time
  s = stks->get(index);

  // loop over all words ending at this node
  if(Debug_flag.extend)
    printf("%d endwords\n", endword->size);
  for(INT j = 0; j < endword->size; j++){
    w = endword->get(j)->id;
    wp = endword->get(j)->prior;
    v = endword->get(j)->version;

    // Don't consider Sentence_end if there's more to come
    if(!Cross_sentence && w == Sentence_end_index && t < proc_frames_seen()) 
      continue;

    //    if(w == Sentence_start_index){
    //      if(!Cross_sentence)
    //	continue;
    //      if(t == sentence_length)
    //	continue;
    //    }
    

    nadded = 0;
    if(hyp_active) assert(hyp_active->size == hyps->size);

    // loop over all hypothesis being propagated
    INT maxoff = INT_MIN;

    // average lm contexts
    INT lmp_init = 0;
    if(average_lm_contexts > 0){
      lmp_init = 0;
      hyps->init_cursor();
      INT hcount = 0;
      while(hyps->next()) {
	hcount++;
	lmp_init += (lm->get_prob(hyps->read(), w, hcount))/(hcount+1);
      }
      //      lmp /= hyps->size;
      lmp_init /= 3;
    }
      

    // Update LUB first
    hyps->init_cursor();
    if(hyp_active) hyp_active->init_cursor();
    INT hcount = -1;
    offp_queue->fast_clear();
    while(hyps->next()){
      hcount++;
      if(hyp_active) hyp_active->next();
      if(hyp_active == NULL || hyp_active->read()) {
	h = hyps->read();
	if(Debug_flag.extend){
	  printf("Extending %d by %d\n", hcount, w&WORD_MASK);
	  h->log_display();
	}
	lmp = lm->get_prob(h, w, hcount);
	if(average_lm_contexts > 0) 
	  lmp += lmp_init/average_lm_contexts;
	offp = wp + lmp - h->get_delta();
	if(Debug_flag.extend)
	  printf(" offp(%d) = wp(%d) + lmp(%d) + delta(%d)\n", 
		 offp, wp, lmp, h->get_delta());
	offp_queue->enqueue(offp);
	maxoff = MAX(maxoff, offp); // update lub at end
      } 
    }
    // This will be an over-estimate for LUB in the case where maxoff
    //  was at a cross-sentence
    if(strict_lub && maxoff != INT_MIN && 
       delta(Bestprob_this_ext+maxoff, Best_t)<0)
      backtrace_lub(t, maxoff);
  
    // Now extend hyps
    hyps->init_cursor();
    if(hyp_active) hyp_active->init_cursor();
    hcount = -1;
    while(hyps->next()){
      hcount++;
      if(hyp_active) hyp_active->next();
      if(hyp_active == NULL || hyp_active->read()) {
	offp = offp_queue->dequeue();
	p = logprob[EXIT_STATE] + offp;

	// t-1 since this frame's acoustic data hasn't been used for exit state
	if(delta(p, t-1) > beam)
	  continue;
	
	h = hyps->read();
	ap = logprob[EXIT_STATE] - h->get_logprob() - h->get_delta();

	if(Debug_flag.extend)
	  printf(" %d:add_hyp(%d, %d, %d = %d + %d)\n", hcount, t, w, p, 
		 logprob[EXIT_STATE], p-logprob[EXIT_STATE]);
	newh = s->add_hyp(h, t, w, p, ap, v, wp, lm);
	if (newh != NULL){
	  nadded++;
	  if(State_decode || Pause_decode)
	    node_sequence(t, offp, newh);
	}
      }
    }
    res += nadded;

  }
  return(res);    
}

VOID
NODE::node_sequence(INT t, INT offp, HYP* h) {
  INT oldt;
  INT endpause = -1;

  NODE *nd = this;
  INT pausep = 0;
  while(nd != NULL) {
    if(t == nd->entry_time){
      nd = nd->get_parent();
      continue;
    }

    if(State_decode) {
      h->get_nodeseq()->push(nd->base_id);
      h->get_nodetimeseq()->push(t);
      h->get_nodeprobseq()->push(nd->exitprobs->get(t - nd->entry_time) + offp);
    }

    // This assumes that previous() hasn't been deleted by online output
    //  - which shouldn't happen as this func is called where h is a newly 
    //   extended hyp
    if(Pause_decode && nd->get_parent() == NULL) {
      assert(h->get_previous() != NULL);
      endpause = t;
      pausep = nd->exitprobs->get(t - nd->entry_time) - 
	        h->get_previous()->get_delta();
      h->insert_pause(endpause, pausep);
    }

    oldt = t;
    t = nd->parent_exit_times->get(oldt - nd->entry_time);

    nd = nd->get_parent();
  }

}

VOID
NODE::backtrace_lub(INT t, INT offp) {
  INT p, oldt;

  NODE *nd = this;
  while(nd != NULL) {
    if(t == nd->entry_time){
      nd = nd->get_parent();
      continue;
    }

    oldt = t;
    t = nd->parent_exit_times->get(oldt - nd->entry_time);

    INT i, indx;
    BACKTRACE_LIST* lptr = nd->maxprobs;
    for(i = t, indx = t-nd->entry_time; i < oldt; i++, indx++) {
      p = lptr->get(indx) + offp;
      if(delta(p,i) < 0)
	update_lub(p, i);
    }

    nd = nd->get_parent();
  }

}

// this is called as pause_model->hmm_set_lub(0) at the beginning of a sentence
VOID
HMM::hmm_set_lub(INT t) {
  BOOL bt = Backtrace;
  Backtrace = FALSE;
  clear_logprob();
  set_init_logprob(0, t);
  while(!is_floor(output_prob(t))) {
    if(step_forward(t) == FALSE)
      break;
    if(logprob[FIRST_STATE] > lub->get(t))
      update_lub(logprob[FIRST_STATE], t);
    t++;
  }
  Backtrace = bt;
}

// returns TRUE if not end-of-sentence
BOOL
HMM::step_forward(INT t){
  INT j, k, p;
  INT *tmp_ptr;;
  HMM_TRANS* tptr;
  INT *lptr, *nlptr;

  // Prop forward from non-exit states (ie all listed transitions)
  memset(newlogprob, 0x7f, nstates*sizeof(INT));
  if(Backtrace) {
    memset(newstatebp, 0x7f, nstates*sizeof(INT));
    for(tptr = trans; tptr < trans+ntrans; tptr++){
      j = tptr->from;
      if(logprob[j] == PROB_NULL)
	continue;
      p = logprob[j] + tptr->prob;
      k = tptr->to;

      nlptr = newlogprob+k;
      if(*nlptr == PROB_NULL) {
	*nlptr = p;
	newstatebp[k] = statebp[j];
      } else if (Forward_process) {
	*nlptr = logadd(p, *nlptr);
	if (*nlptr < p)
	  newstatebp[k] = statebp[j];
      } else if (*nlptr < p) {
	*nlptr = p;
	newstatebp[k] = statebp[j];
      } else {
	continue;
      }
    }
    tmp_ptr = statebp;
    statebp = newstatebp;
    newstatebp = tmp_ptr;
  } else {
    for(tptr = trans; tptr < trans+ntrans; tptr++){
      //    for(i = 0; i < ntrans; i++){
      j = tptr->from;
      if(logprob[j] == PROB_NULL)
	continue;
      p = logprob[j] + tptr->prob;
      k = tptr->to;

      nlptr = newlogprob+k;
      if(*nlptr == PROB_NULL) {
	*nlptr = p;
      } else if (Forward_process) {
	*nlptr = logadd(p, *nlptr);
      } else if (*nlptr < p) {
	*nlptr = p;
      } else {
	continue;
      }
    }
  }
  
  tmp_ptr = logprob;
  logprob = newlogprob;
  newlogprob = tmp_ptr;


  if(Backtrace) {
    if(reached_exit()) 
      parent_exit_times->push(statebp[EXIT_STATE]);
    else
      parent_exit_times->push(-2);
    //    assert(t-entry_time == parent_exit_times->size-1);
    if(Pause_decode || State_decode) 
      exitprobs->push(logprob[EXIT_STATE]);
  }


  p = output_prob(t);

  if(t == get_sentence_length()) 
    return(FALSE);

  // exit state doesn't use acoustics
  maxprob = INT_MIN;
  for(lptr = logprob+FIRST_STATE; lptr < logprob+nstates; lptr++) {
    if (*lptr != PROB_NULL) {
      *lptr +=  p;
      if(Backtrace) 
	maxprob = MAX(*lptr, maxprob);
    }
  }

  if(Backtrace) 
    maxprobs->push(maxprob);
  

  return(TRUE);
}


VOID 
NODE::depthfirst_word(INT t, LIST<HYP*>* hyps, GRAMMAR* lm) {
  NODE *nd = this;
  NODE *nextnd;
  HYP* h;

  while(t < sentence_length) {
    if(nd->step_forward(t) == FALSE)
      break;

    INT maxp = nd->get_maxprob();
    if(nd->reached_exit()) {
      // backtrace from word ends using most probable LM context
      if(t > 0) {
	INT maxoff = INT_MIN;
	for(INT j = 0; j < nd->endword->size; j++) {
	  //	  INT w = nd->endword->get(j);
	  //	  INT wp = nd->wordprior->get(j);
	  INT w = nd->endword->get(j)->id;
	  INT wp = nd->endword->get(j)->prior;
	  INT k = 0;
	  hyps->init_cursor();
	  while(hyps->next()) {
	    h = hyps->read();
	    maxoff = MAX(maxoff, 
			 (wp + lm->get_prob(h, w, k) - h->get_delta()));
	    k++;
	  }
	}
	if(maxoff > INT_MIN) 
	  nd->backtrace_lub(t, maxoff);
      }
      
      // find most probable successor
      INT addp = nd->get_exit_prob() + logcode(0.5);
      nextnd = nd;
      for(INT i = 0; i < nd->succ->size; i++) {
	INT op = nd->succ->get(i)->output_prob(t) + addp;
	if(op > maxp) {
	  maxp = op;
	  nextnd = nd->succ->get(i);
	}
      }

      // if out of beam quit else continue
      if(nd->delta(maxp, t) >= state_beam) {
	break;
      } else if(nd == nextnd) {
	t++;
      } else {
	nd->deactivate();
	nextnd->activate(t);
	nd = nextnd;
      }	
    } else {
      // if out of beam quit else continue
      if(nd->delta(maxp, t) >= state_beam) 
	break;
      else
	t++;
    }
  }
  nd->deactivate();
  return;
}

void* 
BACKTRACE_LIST::operator new(size_t size) {
  if(size != sizeof(BACKTRACE_LIST))
    return ::new char[size];

  BACKTRACE_LIST* res = free_list;

  if(res == NULL){
    // add more memory to the free list

    // use ::new to avoid calling the constructor
    BACKTRACE_LIST* new_block = 
      (BACKTRACE_LIST*) ::new char[block_size * sizeof(BACKTRACE_LIST)];

    if(new_block == NULL) return NULL;
    
    for(INT i = 0; i < block_size-1; i++) {
      new_block[i].previous = &new_block[i+1];
      new_block[i].alloc = DEF_BACKTRACE_LIST_ALLOC;
      new_block[i].default_alloc = DEF_BACKTRACE_LIST_ALLOC;
      new_block[i].data = new INT[DEF_BACKTRACE_LIST_ALLOC];
      new_block[i].size = 0;
    }

    new_block[block_size-1].previous = NULL;
    new_block[block_size-1].alloc = DEF_BACKTRACE_LIST_ALLOC;
    new_block[block_size-1].default_alloc = DEF_BACKTRACE_LIST_ALLOC;
    new_block[block_size-1].data = new INT[DEF_BACKTRACE_LIST_ALLOC];
    new_block[block_size-1].size = 0;
    res = new_block;

    total_backtrace_lists += block_size;

  }
   
  free_list = res->previous;
  return res;  
}

void
BACKTRACE_LIST::operator delete(void* dead, size_t size){
  if(size != sizeof(BACKTRACE_LIST)){
    ::delete [] ((char*) dead);
    return;
  } 
  ((BACKTRACE_LIST*)dead)->previous = free_list;
  ((BACKTRACE_LIST*)dead)->size = 0;
  free_list = (BACKTRACE_LIST*)dead;
}


VOID
BACKTRACE_LIST::free() {
  previous = free_list;
  size = 0;
  free_list = this;
}
