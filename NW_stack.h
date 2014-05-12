// -*- Mode: C++;  -*-
// File: stack.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Stack for stack decoder, derived from priority
//*             queue of hypotheses
//*
//* CLASSES:   STACK
//* 
//* REQUIRES:  PQUEUE (pqueue.h) 
//*            HYP (NW_hypothesis.h)
//*
//* HISTORY:
//*  May 24 15:42 1994 (sjr): Tested okay within decoder.  Beware copy,
//*                            which is not deep...
//* Created: Wed May 18 09:49:47 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



class STACK : public MINPQUEUE<HYP*> {
protected:
  static INT max_hyps;
  static BOOL merge_mode;

 public:
  STACK(INT n = 15) : MINPQUEUE<HYP*>(n) {   }

  // works around compiler errors in 2.6.3 and 2.7.0
  ~STACK() { 
    if(data != NULL) {
      delete[] data;
      data = NULL;
    }  
  }

  INT get_max_hyps() { return(max_hyps); }
  VOID set_max_hyps(INT n) { max_hyps = n; }
  VOID set_merge_mode() { merge_mode = TRUE; }

  VOID merge(INT i, INT p) {
    get(i)->add_to_logprob(p);
    // min-heap property is okay above i, so need to heapify_down
    // to the bottom
    heapify_down(i+1,size);
  }

  //  HYP* add_hyp(HYP *h, LATTICE_LINK *le, INT p, GRAMMAR *lm) {
  //      return(add_hyp(h, le->get_end_time(), le->get_word(), p,  
  //  		   le->get_version(), lm));
  //    }

  HYP* add_hyp(HYP *h, INT t, INT w, INT p, INT wdp, INT v, INT r, 
		      GRAMMAR *lm) {
    INT i;
    HYP *reth = NULL;
    INT ow = w;

    if(lm->use_unk || lm->tagged) 
      ow = lm->oov2unk(w);

    if(Debug_flag.hyp){
      fprintf(stderr, "t = %d, w = %d ow = %d\n", t, w, ow);
      h->info(NULL, NULL, gstderr, "h being extend by STACK::add_hyp()");
    }

    // Don't add it to a full stack where all hypotheses are more
    //  probable
    // Can merge hyps if desired
    if(size >= max_hyps && p <= top()->get_logprob()) {
      return(NULL);
    }

    // Always add to an empty stack
    if(size == 0) {
      if(Debug_flag.hyp)
	fprintf(stderr, "Adding to empty stack\n");
      reth = new HYP(w, t, p, wdp, h, ow, v, r);
      insert(reth);
      return(reth);
    }

    // Add to the stack if newh is either novel hypothesis or has higher 
    //  probability than existing hypothesis

    // If merging check to see if we should perform a merge with a 
    //  hyp in the stack
    if(merge_mode) {
      for(i = size-1; i >= 0; i--){
	// use hyp_exact_similar for merging
	if(lm->hyp_exact_similar(get(i), h, w)){
	  merge(i, p);
	  //	  if(ha != NULL) delete h;
	  return(NULL);
	}
      }
    } else {
      // if not merging then replace an lm-similar hypothesis with lower prob
      for(i = size-1; i >= 0; i--){
	// if we are going to place a new hyp on the stack we
	//  need to check for any hypotheses already on the stack that
	//  are "language model similar"
	if(lm->hyp_similar(get(i), h, w)){
	  if(p > get(i)->get_logprob()) {
	    reth = new HYP(w, t, p, wdp, h, ow, v, r);
	    Hyp_garbage_list->push(swap(i, reth));
	    return(reth);
	  } else {
	    //	    if(ha != NULL) delete h;
	    return(NULL);
	  }
	}
      }
    }

    // h is now a different hypothesis to any on the stack and is if the 
    //  stack is full, h has a better score than the lowest in the stack
    if(Debug_flag.hyp)
      fprintf(stderr, "adding to stack\n");
    reth = new HYP(w, t, p, wdp, h, ow, v, r);
    if(size >= max_hyps) 
      Hyp_garbage_list->push(pop());
    insert(reth);
    return(reth);
  }

};


