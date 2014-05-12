// -*- Mode: C++;  -*-
// File: node.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Node in a lexicon tree, principal compute engine
//*             for decoder.  Derived from HMM phone models, contains
//*             cached LM probs, plus lexical information
//*
//* CLASSES: NODE (nested class of LEXICON)
//* 
//* REQUIRES: HMM, TREE_LEXICON, HYP
//*           #included in lexicon.h
//*
//* HISTORY:
//*  Jan 27 12:58 1995 (sjr): Links back to parent, and state path decoding
//*  May 24 15:45 1994 (sjr): Pruning, forward and create_hyps all test
//*                            okay in decoder
//* Created: Thu Apr 14 10:49:33 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GLOBAL_VAR(INT, Bestprob_this_ext);
GLOBAL_VAR(INT, Best_t);

class WORDINF {
public:
  WORDINF(INT i, INT p, INT v) {
    id = i;
    prior = p;
    version = v;
  }

  INT id;
  INT prior;
  INT version;
};

class NODE : public SUBWORD {
  friend class TREE_LEXICON;
private:
  static INT nnodes;		// total num nodes in tree
  static QUEUE<INT>* offp_queue;
protected:
  // Language model stuff for a set of hypotheses

  
  // relating to the structure and pathway of the tree
  NODE *parent;
  LIST<NODE*> *succ;		// size is num successors
  LIST<INT> *pathword;	        // size is num pathwords
  LIST<BOOL> *hyp_active;	// size is num hypotheses -- dynamic
  LIST<INT> *maxlmprob;	        // size is num hypotheses -- dynamic
  INT maxlmub;			// max im maxlmprob
  INT maxuniprob;
  INT base_id;
  INT level;			// level in tree

  // words ending at this node
  LIST<WORDINF*> *endword;
  //  LIST<INT> *endword;		// size is num words ending at this node 
  //  LIST<INT> *wordprior;         // prior on pronunciation corresponds to endword
  //  LIST<INT> *pronversion;       // pronunciation version corresponds to endword

public:

  NODE(HMM *p, NODE *pp, INT nh = 10, INT hp = TRUE) : SUBWORD(p) {
      //      trans = p->get_trans();
      //      ntrans = p->get_ntrans();
      succ = new LIST<NODE*>(1);
      endword = new LIST<WORDINF*> (1);
      if(smear_unigram || pathword_threshold > 0)
	pathword = new LIST<INT>(1);
      else
	pathword = NULL;
      parent = pp;
      if(hp) {
	maxlmprob = new LIST<INT>(nh+1);
	hyp_active = new LIST<BOOL>(nh+1);
      } else {
	maxlmprob = NULL;
	hyp_active = NULL;
      }
      maxlmub = INT_MIN;
      maxuniprob = INT_MIN;
      retrieved_lmprob = FALSE;
      level = 0;
      //      deactivate();
      nnodes++;
  }

  // returns the successor node if ph matches a node in the succ list,
  //  else returns NULL
  NODE* match_succ(HMM *ph);

  // creates a node corresponding to HMM ph and adds to succ list
  NODE* create_successor(HMM* ph, INT nh = 10, BOOL hp = TRUE);

  VOID node_sequence(INT t, INT p, HYP* h);
  VOID backtrace_lub(INT t, INT p);

  NODE* get_parent() { return(parent); }

  BOOL is_less_than(NODE* nd) { return(maxprob < nd->maxprob); }
  BOOL is_greater_than(NODE* nd) { return(level > nd->level); }

  BOOL is_equal(HMM* ph)
  {  return (nstates==ph->get_nstates() && !strcmp(name, ph->name)); }

  INT nsucc() { return(succ->size); }
  INT npathwords() { return(pathword->size); }

  INT get_beam() { return(beam); }
  // get rid of vacant space in the fixed lists
  VOID squeeze() {
    succ->squeeze();
    if(smear_unigram || pathword_threshold > 0)
      pathword->squeeze();
    endword->squeeze();
    //    wordprior->squeeze();
    //    pronversion->squeeze(); 
  }

  static BOOL smear_unigram;
  static BOOL smear_context;
  BOOL retrieved_lmprob;

  BOOL forward(INT t, LIST<HYP*>* hyps, GRAMMAR* lm, 
	       QUEUE<QUEUE<NODE*>*>* act, INT& maxlevel);
  BOOL prune(INT t, LIST<HYP*> *hyps, GRAMMAR *lm);
  INT  create_new_hyps(INT t, LIST<HYP*> *hyps, GRAMMAR *lm, QUEUE<STACK*> *stks);
  VOID depthfirst_word(INT t, LIST<HYP*>* hyps, GRAMMAR* lm); 


  BOOL activate(INT t, INT p = -1) { 
    p = parent->get_exit_prob();
    BOOL res = SUBWORD::activate(t, p);
    nnodes_used++;

    retrieved_lmprob = FALSE;
    maxlmub = parent->maxlmub;
    if(maxlmprob != NULL) {
      maxlmprob->fast_clear();
      if(parent->nsucc() == 1) {
	if((retrieved_lmprob = parent->retrieved_lmprob) == TRUE)
	  maxlmprob->copy(parent->maxlmprob);
      }
    }

    if(hyp_active != NULL) {
      hyp_active->fast_clear();
      hyp_active->copy(parent->hyp_active);
    }

    return(res);
  }

  VOID short_info(INT t = -1) {
    printf("%s/%d(t=%d,%d): maxprob = %d, maxlmub = %d\n  logprobs: ", 
	   name,id, t, start_time, maxprob, maxlmub);
    printf("entry_time = %d, frame_activated = %d, last_activated = %d\n",
	   entry_time, frame_activated, last_activated);
    if(Backtrace) {
      printf("maxprobs[%d]  parent_exit_times[%d]\n",
	     maxprobs->size, parent_exit_times->size);
      if(State_decode || Pause_decode)
	printf("exitprobs[%d]\n", exitprobs->size);
    }
    for(INT i = 0; i < nstates; i++)
      if(logprob[i] == PROB_NULL)
	printf("NULL ");
      else
	printf("%d ", logprob[i]);
    if(t == -1)
      printf("\n");
    else if(t < proc_frames_seen())
	printf("lub[%d] = %d\n", t, lub->get(t));
    fflush(stdout);
  }

  VOID info(INT t = -1, STRING s = "") {
    printf("t = %d (start_time = %d): %s\n(%s,%d) ", 
	   t, start_time, s, name, id);
    if(activated)
      printf("***ACTIVATED***   ");
    else
      printf("++INACTIVE++   ");
    if(smear_unigram || pathword_threshold > 0)
      printf("%d successors, %d pathwords, %d endwords, lub = %d\n", 
	     nsucc(), pathword->size, endword->size, 0);//lub->get(t));
    else
      printf("%d successors, %d endwords, lub = %d\n", 
	     nsucc(), endword->size, 0);//lub->get(t));
    if(logprob[0] == PROB_NULL)
      printf("  logprob = (VOID, ");
    else
      printf("  logprob = (%d/%f, ", logprob[0], logdecodelog(logprob[0]));
    for(INT i = 2; i < get_nstates(); i++)
      if(logprob[i] == PROB_NULL)
	printf("VOID, ");
      else
	printf("%d/%.3f, ", logprob[i], logdecodelog(logprob[i]));
    if(logprob[1] == PROB_NULL)
      printf("VOID)\n");
    else
      printf("%d/%.3f)\n", logprob[1], logdecodelog(logprob[1]));
    printf(")\n");
  }

};
