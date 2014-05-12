 // -*- Mode: C++;  -*-
// File: lexicon.cc
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* RELATED PACKAGES:
//*
//* HISTORY:
//*  Feb 21 12:50 1995 (sjr): *_in() calls now include lub_n parm
//*  Apr 12 14:09 1994 (sjr): Tree-structuring in place, tested build and info
//* Created: Fri Apr  8 09:57:31 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "all.h"



VOID 
LINEAR_LEXICON::build(PARAM_TABLE *parm, const INT nw, BOOL phi, GENFILE *lexfp, 
		      GENFILE *phonfp, GENFILE *priorfp) {

  //  BOOL add_sil = FALSE;
  CHAR buf[BUFSIZ];
  CHAR wd[BUFSIZ];
  CHAR lastwd[BUFSIZ];
  CHAR *pptr, *eptr, *pronptr = NULL;
  CHAR *endp = new CHAR;
  FLOAT pronprior = 1.0;
  INT logpronprior = 0;
  //  INT version = 0;
  BOOL novel;
  INT word_id, last_word_id;
  INT wordcount;
  INT i;
  BOOL ppflag = (!parm->get_bool("no_pron_prior"));
  WORD *word_model = NULL;
  
    // build the phone models and store in a hash table
  LIST<FLOAT>* priors = make_phoneset(parm, phi, phonfp, priorfp);
  if(Verbose){
    printf("Built phoneset, %d phones\n", phoneset->size());
    fflush(stdout);
  }

  // set the root node to an optional interword pause
  pause_model = phoneset->get("interword-pause");
  if(pause_model == NULL)
    panic("Phoneset must include interword-pause model");
  Silence_distr = pause_model->get_distr();

  wordcount = npron  = 0;

  if(Silence_index == -1) {
    // Silence not in LM so not yet in vocab
    Silence_index = vocabulary->size;
    vocabulary->push(Silence_symbol);
    vocab_hash->insert(Silence_symbol, Silence_index);
    PRONUNCIATION* p = new PRONUNCIATION();
    p->get_subwords()->enqueue(new SUBWORD(pause_model));
    WORD* w = new WORD(Silence_index);
    w->get_prons()->enqueue(p);
    get_words()->put_and_grow(Silence_index, w, NULL);
    wordcount++;
    npron++;
  }

  BOOL *checkwords = new BOOL[nw]; 
  for(i = 0; i < nw; i++) checkwords[i] = FALSE;
  INT dudphone = 0;
  last_word_id = -1;

  while(lexfp->fgets(buf, BUFSIZ) != NULL) {
    // miss out comments
    if(buf[0] == '#')
      continue;

        // word string
    if((eptr = strpbrk(buf+1, "( \t")) == NULL)
      panic("Failed to read word %d\n%s", npron, buf);
    strncpy(wd, buf, (eptr - buf)); wd[eptr-buf] = '\0';
    assert(strlen(wd) > 0);

    // get word index
    word_id = vocab_hash->get(wd);
    if(word_id == -1 ) 
      continue;

    strcpy(lastwd, wd);

    // prior on pronunciation - ignore for now!

    if((pptr = strchr(buf+1, '(')) != NULL){
      pptr++;
      // pptr points to prior on pronunciation
      pronprior = strtod(pptr, &endp);
      if(endp == pptr)
	panic("Failed to read pronunciation %d\n%s", npron, buf);
      if(*endp == ')')
	pronptr = endp + 1;			// read correctly
      else
	panic("Expected ')' in pronunciation %d\n%s", npron, buf);
    
      // Don't need to include 0 prob pronunciations
      if(ppflag) {
	if(pronprior <= 0.0) {
	  if(Debug_flag.tree) {
	    printf("0 prob pronunciation for %s\n %s\n", wd, buf);
	    fflush(stdout);
	  }
	  continue;
	}
	logpronprior = logcode(pronprior);
      } else {
	logpronprior = 0;
      } 
    } else {
      pronptr = eptr;
      logpronprior = 0;
    }

    // Build the word models
    novel = (word_id != last_word_id);
    if(novel) {
      assert(word_id < nw);
      checkwords[word_id] = TRUE;
      ++wordcount;
      word_model = new WORD(word_id);
      get_words()->put_and_grow(word_id, word_model, NULL);
      last_word_id = word_id;
    }
    ++npron;

    // optional silence at the start of each word
    SUBWORD* next_sw = new SUBWORD(pause_model);
    PRONUNCIATION* pron = new PRONUNCIATION();
    pron->get_subwords()->enqueue(next_sw);
    HMM* next_phone;
    while((pptr = strtok(pronptr, " \t\n")) != NULL) {
      pronptr = NULL;
      if(!strcmp(pptr, "+"))
	continue;
      next_phone = phoneset->get(pptr);
      if(next_phone == NULL){
	dudphone++;
	fprintf(stderr, "LINEAR_LEXICON::build(), illegal phone %s in word %s\n",
	      pptr, wd);
	break;
      }
      next_sw = new SUBWORD(next_phone);
      pron->get_subwords()->enqueue(next_sw);
    }
    word_model->get_prons()->enqueue(pron);
  }
  if(dudphone > 0)
    panic("(At least) %d bogus phones in pronunciations\n", dudphone);
  
  if(wordcount < nw){
    for(i = 0; i < nw; i++)
      if(checkwords[i] == FALSE) 
	fprintf(stderr, "%s\n", vocabulary->get(i));
    fprintf(stderr, 
	    "Not all words in LM are present in dictionary (wordcount = %d)\n", 
	    wordcount);
  }
  delete checkwords;

  set_acoustic_input(parm, priors);

  delete priors;

  if(Verbose || Debug_flag.tree_stats){
    printf("Built dictionary\n");
    //    info();
    printf("beam = %.2f,  state_beam = %.2f, phone threshold = %f\n", 
	   logdecodelog(pause_model->get_beam()), 
	   logdecodelog(pause_model->get_state_beam()), 
	   pause_model->get_lna_floor());
  }
}

VOID 
LINEAR_LEXICON::score(HYP* , GRAMMAR*, LIST<INT>*, LIST<INT>*){
  panic("LINEAR_LEXICON::score() not implemented\n");
}


INT
LINEAR_LEXICON::extend(INT , LIST<HYP*>*, GRAMMAR*, QUEUE<STACK*>*){
  panic("LINEAR_LEXICON::extend() not implemented\n");
  return 0; 
}


// Arguments  t: new word start time
//        hlist: list of HYPs to extend
//           lm: GRAMMAR probably finite state network
//         stks: QUEUE  of stacks to insert word-extended HYPs
// Operation:   Extend all HYPs in hlist by all word extensions permitted by lm
//              inserting extend HYPS into stks
// Return:   number of extended HYPs created
INT 
LINEAR_LEXICON::align(INT t, LIST<HYP*> *hlist, GRAMMAR *lm, QUEUE<STACK*> *stks){
  LIST<INT>* slist;
  HYP* h;
  WORD* w;
  INT res = 0;
  hlist->init_cursor();
  while(hlist->next()) {
    h = hlist->read();
    slist = lm->get_successors(h);
    slist->init_cursor();
    while(slist->next()){
      w = words->get(slist->read());
      res += w->extend(t, h, stks);
    }
  }
  return res;
}

// Arguments  t: word start time
//            h: HYP to extend
//         stks: QUEUE  of stacks to insert word-extended HYPs
// Operation:   Extend h by by this WORD, inserting extended HYPS in stks
// Return:   number of extended HYPs created
INT
WORD::extend(INT startt, HYP* h, QUEUE<STACK*> *stks){
  INT res = 0;
  PRONUNCIATION* pr;
  SUBWORD* sw;
  INT v = 0;

  get_prons()->init_cursor();
  while(get_prons()->next()){
    pr = get_prons()->read();

    // activate the first phone
    pr->get_subwords()->top()->activate(startt, h->get_logprob());
    for(INT t = startt; t <= pr->get_subwords()->top()->get_sentence_length(); t++) {
      pr->get_subwords()->init_cursor();
      while(pr->get_subwords()->next()) {
	sw = pr->get_subwords()->read();
	if(!sw->forward(t))
	  sw->deactivate();	
      }
      sw = pr->get_subwords()->last();
      if(sw->logprob[EXIT_STATE] != PROB_NULL && t > 0){
	//create_new_hyp(t, h, stks);
	INT index = t - (startt + 1);
	if(index < 0) {
	  if(t == startt)
	    return(0);
	  else
	    panic("t = %d, start_time = %d: Rogue hypothesis in word %d\n", 
		  t, startt, word_id);
	}

	// create enough new stacks to reach the ref time of the new hypotheses
	while(index >= stks->size)
	  stks->enqueue(new STACK());
	
	// get the stack with the relevant ref time
	INT prob = h->get_delta() + sw->logprob[EXIT_STATE];
	if(sw->delta(prob, t-1) < sw->get_beam())
	  stks->get(index)->insert(new HYP(get_word_id(), t, prob, prob-h->get_logprob(),  h, get_word_id(), v));
      }
    }

    // Finally make sure everything is deactivated
    pr->get_subwords()->init_cursor();
    while(pr->get_subwords()->next()) 
      pr->get_subwords()->read()->deactivate();
    v++;
  }

  return res;
}
