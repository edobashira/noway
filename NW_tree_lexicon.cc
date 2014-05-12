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
TREE_LEXICON::init_extend(INT t, LIST<HYP*> *hyps, GRAMMAR *) {
  HYP *h;
  INT p;

  // root has no parent node, so we can't call root->activate()
  root->clear_logprob();
  root->set_init_logprob(hyps->top()->get_logprob(), t);
  root->set_startprob(hyps->top()->get_logprob());
  root->activated = TRUE;
  root->nnodes_used++;
  // set LM probs to default
  root->retrieved_lmprob = FALSE;
  root->maxlmub = INT_MIN;

  root->entry_time = t;
  root->frame_activated = t;
  root->last_activated = root->sentence_count;
  if(Backtrace) {
    root->parent_exit_times->fast_clear();
    root->maxprobs->fast_clear();
    if(State_decode || Pause_decode)
      root->exitprobs->fast_clear();
    assert(root->get_active_sw_list()->size == 0);
  }

  if(root->hyp_active != NULL) {
    root->hyp_active->fast_clear();
    root->hyp_active->set(TRUE, hyps->size);
  }

  if(root->maxlmprob != NULL) {
    root->maxlmprob->fast_clear();
    hyps->init_cursor();
    while(hyps->next()){
      h = hyps->read();
      //    p = lm->get_default_prob(h->history()) - h->get_delta();
      p = -(h->get_delta());
      root->maxlmprob->push(p);
    }
    if(root->hyp_active != NULL) {
      assert(root->maxlmprob->size == root->hyp_active->size);
      assert(root->hyp_active->size == hyps->size);
    }
  }
  root->maxlmub = 0;
}

// Given a list of hypotheses starting at time t, this function
//  performs a viterbi through the tree, adding new hypotheses to
//  the appropriate stacks.  Returns number of new hypotheses created.
INT
TREE_LEXICON::extend(INT t, LIST<HYP*> *hyps, GRAMMAR *lm, 
		  QUEUE<STACK*> *stks) {
  NODE* nd;
  QUEUE<QUEUE<NODE*>*>* tmp_active;
  INT res = 0;
  INT i;
  BOOL activedirty = FALSE;

//   if(root->strict_lub) {
//    // save state
//    INT sc = root->sentence_count;
//    BOOL sd  = State_decode;
//    BOOL pd = Pause_decode;
//    root->sentence_count = -sc;
//    State_decode = FALSE;
//    Pause_decode = FALSE;

//    init_extend(t, hyps, lm);
//    root->depthfirst_word(t, hyps, lm);

//    while(root->get_active_sw_list()->size > 0)
//      root->get_active_sw_list()->pop()->backtrace_garbage_collect();

//    // restore state
//    State_decode = sd;
//    Pause_decode = pd;
//    root->sentence_count = sc;
//  }

  init_extend(t, hyps, lm);

  // Carry on extending hypotheses while nodes are active and the sentence
  //  hasn't finished
  //  active->fast_clear();
  //  active->push(root); 
  //  next_active->fast_clear();


  //  for(i = 0; i < active->size; i++)
  //    assert(active->get(i)->size == 0);

  //  for(i = 0; i < next_active->size; i++)
  //    assert(next_active->get(i)->size == 0);

  active->get(0)->push(root);

  INT maxlevel = 1;

  if(Debug_flag.tree_stats) {
    for(i = 0; i < depth_count->size; i++)
      depth_count->put(i, 0);
    pwcount = 0;
  }

  // t is an arg to this function, node_count is just for debug purposes
  for(node_count = 0; t <= sentence_length(); t++) {
    Bestprob_this_ext = INT_MIN;
    // loop over active nodes at time t
    INT next_active_size = 0;
    for(INT level = 0; level <= maxlevel; level++) {
      QUEUE<NODE*>* nds = active->get(level);
      while(nds->size > 0) {
	nd = nds->pop();
	node_count++;
	if(Debug_flag.tree_stats) {
	  depth_count->put(level, (depth_count->get(level) + 1));
	  if(nd->npathwords() > nd->pathword_threshold)
	    pwcount++;
	}
	if(nd->forward(t, hyps, lm, active, maxlevel)) {
	  next_active->get(level)->push(nd);
	  next_active_size++;
	} else {
	  nd->deactivate();
	}

	// need to check since this may be the first time we hit end-of-sentence
	if(t > sentence_length()){
	  activedirty = TRUE;
	  break;
	}

	if(nd->endword->size > 0 && nd->logprob[EXIT_STATE] != PROB_NULL && t > 0)
	  res += nd->create_new_hyps(t, hyps, lm, stks);
      }
    }
    // need to check since this may be the first time we hit end-of-sentence
    if(t > sentence_length())
      break;
    
    // active <-- next_active
    tmp_active = active;
    active = next_active;
    next_active = tmp_active;  // should be clear

    if(next_active_size == 0)
      break;
  }

  if(activedirty)
    for(i = 0; i < active->size; i++)
      while(active->get(i)->size > 0)
	active->get(i)->pop()->deactivate();

  if(Backtrace){
    while(root->get_active_sw_list()->size > 0)
      root->get_active_sw_list()->pop()->backtrace_garbage_collect();
  }

  if(Debug_flag.decode_stats)
    Decode_inf_file->fprintf("%d %d\n", res, node_count);
  return(res);
}


// Use the stdio BUFSIZ (1024)
// const INT BUFSIZ = 256;
// Lexicon file format
// <word-string> <optional prior prob> <phone sequence>
VOID
TREE_LEXICON::build(PARAM_TABLE *parm, const INT nw, BOOL phi,
	       GENFILE *lexfp, GENFILE *phonfp, GENFILE *priorfp){
  // one pronunciation per line
  CHAR buf[BUFSIZ];
  CHAR wd[BUFSIZ];
  CHAR lastwd[BUFSIZ];
  CHAR *pptr, *eptr, *pronptr = NULL;
  CHAR *endp = new CHAR;
  FLOAT pronprior = 1.0;
  INT logpronprior = 0;
  INT version = 0;
  BOOL novel;
  INT word_id, last_word_id;
  BOOL new_word;
  INT wordcount, totalwordcount;
  INT nh = parm->get_int("n_hyps", 15);
  BOOL hp = parm->get_bool("hyp_prune");
  BOOL ppflag = (!parm->get_bool("no_pron_prior"));
  BOOL use_unk = parm->get_bool("use_unk") && !parm->get_bool("tagprobs");
  INT i;
  LIST<FLOAT> *priors;
  BOOL pt = FALSE;
  WORDINF *wi;
    

  if(hp && NODE::pathword_threshold == 0) {
    NODE::pathword_threshold = 1;
  }
  
  if(parm->get_bool("tagprobs")) {
    NODE::smear_unigram = FALSE;
    NODE::smear_context = FALSE;
  } else {
    NODE::smear_unigram = parm->get_bool("smear_unigram") || 
	                  parm->get_bool("demo");
    if(Verbose && NODE::smear_unigram && parm->get_bool("demo"))
      printf("-demo option setting -smear_unigram\n");
    NODE::smear_context = parm->get_bool("smear_context");
  }

  if(NODE::smear_unigram || NODE::pathword_threshold > 0) pt = TRUE;

  // build the phone models and store in a hash table
  priors = make_phoneset(parm, phi, phonfp, priorfp);
  if(Verbose){
    printf("Built phoneset, %d phones\n", phoneset->size());
    fflush(stdout);
  }

  // set the root node to an optional interword pause
  pause_model = phoneset->get("interword-pause");
  if(pause_model == NULL)
    panic("Phoneset must include interword-pause model");
  Silence_distr = pause_model->get_distr();

  root = new NODE(pause_model, NULL, nh, hp);
  if(Backtrace) {
    root->parent_exit_times = new BACKTRACE_LIST();
    root->maxprobs = new BACKTRACE_LIST();
    if(State_decode || Pause_decode)
      root->exitprobs = new BACKTRACE_LIST();
  }
  root->base_id = pause_model->get_id();

  nlevels = 1;
  active->enqueue(new QUEUE<NODE*>());
  next_active->enqueue(new QUEUE<NODE*>());

  if(Silence_index == -1) {
    // Silence not in LM so not yet in vocab
    Silence_index = vocabulary->size | LM_OOV_MASK;
    vocabulary->push(Silence_symbol);
    vocab_hash->insert(Silence_symbol, Silence_index);
  }

  INT minuni = 1.0;
  if(NODE::smear_unigram) {
    for(i = 0; i < unigrams->size; i++){
      root->maxuniprob = MAX(root->maxuniprob, unigrams->get(i));
      if(unigrams->get(i) < minuni)
	minuni = unigrams->get(i);
    }
  }


  wordcount = npron = totalwordcount = 0;

  // Sentence only has acoustic realisation in x-sentence case
  wordcount++;
  
  // Sentence end has an acoustic realisation as interword-pause
  if(pt)
    root->pathword->push(Sentence_end_index);
  wi = new WORDINF(Sentence_end_index, 0, 0);
  root->endword->push(wi);

  npron++;
  totalwordcount++;
  wordcount++;
  last_word_id = Sentence_end_index;
  strcpy(lastwd, Sentence_end_symbol);
  
  HMM *lpse = phoneset->get("long-pause");
  // use Sentence_start_index for x-word
  if(lpse != NULL && Cross_sentence) {
    if(Verbose)
      printf("Using long-pause model\n");
    NODE *gap_node = root->create_successor(lpse, nh, hp);
    if(pt)
      gap_node->pathword->push(Sentence_start_index);
    wi = new WORDINF(Sentence_start_index, 0, 0);
    gap_node->endword->push(wi);
    npron++;
    totalwordcount++;
    wordcount++;
    last_word_id = Sentence_start_index;
    strcpy(lastwd, Sentence_start_symbol);
  }

  BOOL *checkwords = new BOOL[nw]; 
  for(i = 0; i < nw; i++) checkwords[i] = FALSE;
  checkwords[Sentence_start_index] = TRUE;
  checkwords[Sentence_end_index] = TRUE; 
  INT dudphone = 0;

  if(Debug_flag.tree) {
    printf("Tree debugging on (nw = %d, wordcount = %d, npron = %d)\n",
	   nw, wordcount, npron);
    fflush(stdout);
  }

  //  while(fgets(buf, BUFSIZ, lexfp) != NULL && wordcount < nw) {
  while(lexfp->fgets(buf, BUFSIZ) != NULL) {
    // miss out comments
    if(buf[0] == '#')
      continue;

    // word string
    // special case for things like "(PAREN"
    if(buf[0] == '(') 
      eptr = strpbrk(buf+1, "( \t");
    else
      eptr = strpbrk(buf, "( \t");

    if(eptr == NULL)
      panic("Failed to read word %d\n%s", npron, buf);
    strncpy(wd, buf, (eptr - buf)); 
    wd[eptr-buf] = '\0';
    assert(strlen(wd) > 0);

    // get word index
    word_id = vocab_hash->get(wd);
    if(word_id == -1) {
      // word not in LM
      // For large vocab this assumes that multiple pronunciations
      //  are contiguous in the dictionary file
      if(use_unk) {
	word_id = vocabulary->size | LM_OOV_MASK;
	vocabulary->push(strdup(wd));
	vocab_hash->insert(strdup(wd), word_id);
	new_word = TRUE;
      } else {
	continue;
      }
    } else {
      new_word = FALSE;
    }

    // This will happen if a word has been added to the vocab by something
    // other than the LM (and -use_unk is not specified).  Currently, this
    // only happens if Silence_symbol is in the dictionary and not the LM...
    if (!use_unk && (word_id&WORD_MASK) >= nw)
      continue;
      
  

    strcpy(lastwd, wd);

    if(word_id == Sentence_start_index || word_id == Sentence_end_index)
      continue;			// Dealt with these already

    // prior on pronunciation

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
    if(logpronprior != 0)
      Use_pron_priors = TRUE;
  
    // Store the pronunciation in the tree creating new nodes if necessary
    novel = (word_id != last_word_id);
    if(novel){
      version = 0;
      if(pt)
	root->pathword->push(word_id);
      if(new_word == FALSE) {
	wordcount++;
	assert((word_id&WORD_MASK) < nw);
	assert(checkwords[word_id&WORD_MASK] == FALSE);
	checkwords[word_id&WORD_MASK] = TRUE;
      }
      totalwordcount++;
      last_word_id = word_id;
    } else {
      version++;
    }
    ++npron;

    NODE* current_node = root;
    INT level = 0;
    NODE* next_node;
    HMM* next_phone;
    while((pptr = strtok(pronptr, " \t\n")) != NULL) {
      pronptr = NULL;
      if(!strcmp(pptr, "+"))
	continue;
      next_phone = phoneset->get(pptr);
      if(next_phone == NULL){
	dudphone++;
	fprintf(stderr, "TREE_LEXICON::build(), illegal phone %s in word %s\n",
	      pptr, wd);
	break;
      }
      level++;
      if(level >= nlevels) {
	nlevels++;
	active->enqueue(new QUEUE<NODE*>());
	next_active->enqueue(new QUEUE<NODE*>());
      }
      if((next_node = current_node->match_succ(next_phone)) == NULL) {
	next_node = current_node->create_successor(next_phone, nh, hp);
	next_node->base_id = next_phone->get_id();
	next_node->level = level;
      } else {
	assert(level == next_node->level);
      }
      current_node = next_node;
      nlinear_nodes++;
      if(pt && (novel || current_node->pathword->size == 0 || 
	 current_node->pathword->top() != word_id))
	current_node->pathword->push(word_id);
      if(NODE::smear_unigram) {
	if (word_id & LM_OOV_MASK)
	  //	  current_node->maxuniprob = MAX(current_node->maxuniprob, unigrams->get(Unknown_word_index));
	  current_node->maxuniprob = MAX(current_node->maxuniprob, minuni);
	else
	  current_node->maxuniprob = MAX(current_node->maxuniprob, unigrams->get(word_id));
      }
    }

    wi = new WORDINF(word_id, logpronprior, version);
    current_node->endword->push(wi);
    //    current_node->endword->push(word_id);
    //    current_node->wordprior->push(logpronprior);
    //    current_node->pronversion->push(version); 

    if(Debug_flag.tree && !(totalwordcount%50000)) {
      printf("Read in %d words (last was %s)\n", totalwordcount, wd);
      fflush(stdout);
    }
  }
  if(Debug_flag.tree && !lexfp->eof()) {
    printf("Remaining in lexicon file:\n");
    while(lexfp->fgets(buf, BUFSIZ) != NULL)
      printf("  %s\n", buf);
    fflush(stdout);      
  }

  if(dudphone > 0)
    panic("(At least) %d bogus phones in pronunciations\n", dudphone);

  // It doesn't matter if some words in the LM do not have pronunciations -
  //  they will just never be recognized
  if(wordcount < nw-1 ||
     (use_unk == FALSE && wordcount < nw) ||
     (use_unk == TRUE && wordcount < nw-1 && 
      checkwords[Unknown_word_index] == FALSE)) {
    nmissing = nw-wordcount;
    GENFILE* mpf;
    if(parm->get_bool("missing_pron_file"))
      mpf = new GENFILE(parm->get_string("missing_pron_file"), "w");
    else
      mpf = gstderr;
    for(i = 0; i < nw; i++) {
      if(checkwords[i] == FALSE) 
	mpf->fprintf("%s\n", vocabulary->get(i));
    }
    delete mpf;
    fprintf(stderr, "Missing LM words in dictionary (wordcount=%d, %d missing pronunciations)\n", wordcount, nmissing);
  } else {
    nmissing = 0;
  }
		 
  delete checkwords;
  
  // maxlmprob
  if(!hp && root->pathword_threshold > 0)
    allocate_maxlmprob(nh);

  // Set acoustic input
  if(!parm->get_bool("no_decode"))
    set_acoustic_input(parm, priors);

  delete priors;

  if(Verbose || Debug_flag.tree_stats){
    printf("Built dictionary\n");
    info();
    printf("beam = %.2f,  state_beam = %.2f, phone threshold = %f\n", 
	   logdecodelog(root->beam), logdecodelog(root->state_beam), 
	   root->lna_floor);
  }

  INT skip_sent = parm->get_int("skip_sents", 0);
  if(Verbose && skip_sent > 0)
    printf("Skipping first %d utterances\n", skip_sent);
  while(skip_sent > 0) {
    root->skip_sentence_acoustics();
    skip_sent--;
  }


  INT skip_frame = parm->get_int("skip_frames", 0);
  if(skip_frame > 0) {
    if(Verbose)
      printf("Skipping first %d frames in first utterance\n", skip_frame);
    if(root->skip_frames_in_current_sent(skip_frame) > 0)
      fprintf(stderr, "Fewer than %d frames in first sentence, so skipped entire sentence\n", skip_frame);
  }

  
  // This just causes memory fragmentation, slowing things down and 
  // increasing process size!
//  commission();

}



// Squeeze list and set up language model lists
VOID
TREE_LEXICON::commission() {
  INT i;
  QUEUE<NODE*> *nodeq = new QUEUE<NODE*>();
  NODE* nd;

  nodeq->push(root);
  while(!nodeq->is_empty()){
    // get the current node
    nd = nodeq->pop();

    // add successor nodes to the queue
    for(i = 0; i < nd->nsucc(); i++)
      nodeq->push(nd->succ->get(i));

    // Delete extraneous space in lists
    nd->squeeze();
  }
  delete nodeq;
}

VOID
TREE_LEXICON::allocate_maxlmprob(INT nh) {
  INT i;
  QUEUE<NODE*> *nodeq = new QUEUE<NODE*>();
  NODE* nd;
  INT nps;

  nodeq->enqueue(root);
  while(!nodeq->is_empty()){
    nd = nodeq->dequeue();
    if(nd->pathword_threshold == 0)
      nps = 0;
    else
      nps = nd->pathword->size;
    if(nps <= nd->pathword_threshold && 
       nd->maxlmprob == NULL)
      nd->maxlmprob = new LIST<INT>(nh);
    for(i = 0; i < nd->nsucc(); i++)
      nodeq->enqueue(nd->succ->get(i));
  }
  delete nodeq;
}

// used in dump lub
VOID
TREE_LEXICON::score(HYP *h, GRAMMAR *lm, LIST<INT> *res, LIST<INT> *res2) {
  INT start, end;
  INT initp;
  HMM *hmm;
  HYP *hh;
  QUEUE<INT> *nts;
  QUEUE<INT> *ns;
  LIST<INT> *endlist = new LIST<INT>();
  LIST<INT> *indxlist = new LIST<INT>();

  // first the acoustic probs
  res->clear();
  BOOL bt = Backtrace;
  Backtrace = FALSE;
  hh = h;
  while(hh != NULL) {
    nts = hh->get_nodetimeseq();
    ns = hh->get_nodeseq();
    assert(nts->size == ns->size);
    for(INT i = 0; i < nts->size; i++) {
      endlist->push(nts->get(i));
      indxlist->push(ns->get(i));
    }
    hh = hh->get_previous();
  } 

  start = 0;    
  initp = 0;
  endlist->init_cursor();
  indxlist->init_cursor();
  while(endlist->next() && indxlist->next()) {
    end = endlist->read();
    hmm = phoneset->get(phonenames->get(indxlist->read()));
    hmm->clear_logprob();
    hmm->set_init_logprob(initp, start);
    for(INT t = start; t <= end; t++) {
      hmm->step_forward(t);
      INT ii = MIN((t-start) + FIRST_STATE, hmm->get_nstates()-1);
      res->push(hmm->get_state_prob(ii));
    }
    initp = hmm->get_exit_prob();
    hmm->clear_logprob();
    res->pop();
    start = end;
  }
  delete endlist;
  delete indxlist;
  Backtrace = bt;

  // Add in the LM probs
  if(lm != NULL) {
    hh = h;
    LIST<INT> *wlist = new LIST<INT>();
    LIST<INT> *tlist = new LIST<INT>();
    while(hh != NULL) {
      wlist->push(hh->get_word());
      tlist->push(hh->get_time());
      hh = hh->get_previous();
    }
    INT w, p;
    INT lmp = 0;
    INT lmp2 = 0;
    LIST<INT> *hist = new LIST<INT>();
    if(res2 != NULL)
      res2->copy(res);
    start = 0;
    wlist->init_cursor();
    tlist->init_cursor();
    while(wlist->next() && tlist->next()){
      w = wlist->read();
      end = tlist->read();
      if(w != Sentence_start_index) {
	if(res2 != NULL)
	  lmp2 = lmp;
	lmp += lm->get_prob(hist, w);
	for(INT k = start; k < end; k++) {
	  p = res->get(k)+lmp;
	  res->put(k, p);
	  if(res2 != NULL) {
	    p = res2->get(k)+lmp2;
	    res2->put(k, p);
	  }
	}
      }
      hist->push(w);
      start = end;
    }
  delete hist;
  delete wlist;
  delete tlist;
  }
  return;
}

VOID 
TREE_LEXICON::info(){
  printf("%d base phones, %d nodes, %d phones in dictionary\n%d words, %d pronunciations, %d missing\n", 
	 phoneset->size(), root->nnodes, nlinear_nodes, vocabulary->size, npron, nmissing);

  if(Debug_flag.tree_stats) {
    INT i;
    QUEUE<NODE*> *nodeq = new QUEUE<NODE*>();
    QUEUE<NODE*> *nextq = new QUEUE<NODE*>();
    QUEUE<NODE*> *tmpq;
    NODE* nd;
    GENFILE* tree_inf_file = new GENFILE("TREE.STATS", "w");
    QUEUE<INT> *pwq = new QUEUE<INT>();

    nodeq->enqueue(root);
    INT depth = 0;
    INT dcount = 0;
    INT pwc = 0;
    INT pwc2 = 0;
    INT nint = 0;
    INT next = 0;
    INT nboth = 0;
    INT pt = NODE::pathword_threshold;
    while(1) {
      dcount = 0;
      pwc2 = 0;
      while(!nodeq->is_empty()){
	nd = nodeq->dequeue();
	if(pt > 0)
	  tree_inf_file->fprintf("%d %d %d %d %d\n", 
		depth, nd->nstates, nd->nsucc(), 
		nd->pathword->size, nd->endword->size);
	else
	  tree_inf_file->fprintf("%d %d %d %d\n", 
		depth, nd->nstates, nd->nsucc(), 
		nd->endword->size);
	  
	if(nd->endword->size == 0)
	  nint++;
	else if(nd->nsucc() == 0)
	  next++;
	else 
	  nboth++;
	dcount++;
	if(nd->npathwords() > pt) {
	  pwc++;
	  pwc2++;
	}
	for(i = 0; i < nd->nsucc(); i++)
	  nextq->enqueue(nd->succ->get(i));
      }
      if(nextq->is_empty()) {
	break;
      }

      // swap nodeq and nextq
      tmpq = nodeq;
      nodeq = nextq;
      nextq = tmpq;
      tmpq = NULL;
      depth_count->enqueue(dcount);
      pwq->enqueue(pwc2);
      depth++;
    }

    delete nodeq;
    delete nextq;
    delete tree_inf_file;
    printf("%d internal nodes, %d external nodes, %d both internal and external\n",
	   nint, next, nboth);
    depth = 0;
    while(depth_count->size > 0) {
      printf("Depth = %d:  %d nodes (%d above pathword threshold)\n", 
	     depth, depth_count->dequeue(), pwq->dequeue()); 
      depth++;
    }
    delete pwq;
    for(i = 0; i < depth; i++)
      depth_count->enqueue(0);    
    printf("%d nodes (%.1f%%) above pathword-threshold (%d)\n", pwc,
	   (FLOAT)pwc/(FLOAT)(nint+next+nboth), pt);
  }
}

// Breadth-first search through the tree
VOID
TREE_LEXICON::dump() {
  INT i;
  QUEUE<NODE*> *nodeq = new QUEUE<NODE*>();
  NODE* nd;
  info();
  nodeq->enqueue(root);
  while(!nodeq->is_empty()){
    nd = nodeq->dequeue();
    for(i = 0; i < nd->nsucc(); i++)
      nodeq->enqueue(nd->succ->get(i));
    printf("Node %d (%s,%d), %d states, %d successors\n", 
	   nd->id, nd->name, nd->distr, nd->nstates, nd->nsucc());
//    printf("( ");
//    for(i = 0; i < nd->pathword->size; i++)
//      printf("%s ", vocabulary->get(nd->pathword->get(i)));
//    printf(")\n");
    if(nd->endword->size != 0){
      printf("    end words: %s (%d), p = %f\n",
	       vocabulary->get(nd->endword->get(0)->id), nd->endword->get(0)->id, 
	       logdecode(nd->endword->get(0)->prior));
      for(i = 1; i < nd->endword->size; i++){
	printf("               %s (%d), p = %f\n",
	       vocabulary->get(nd->endword->get(i)->id), nd->endword->get(i)->id, 
	       logdecode(nd->endword->get(i)->prior));
      }
    }
  }
  delete nodeq;
}

#if 0
// Vile hack required for broken templates in gcc-2.6.3
STR_HASH_TABLE<HMM*>::STR_HASH_TABLE() { 
  tbl = new HASH_TABLE<STR_HASH_ENTRY <HMM*> >(); 
  entry = new STR_HASH_ENTRY <HMM*> ("", (HMM*)NULL);
}

HMM* 
STR_HASH_TABLE<HMM*>::get(STRING s) { 
  entry->key = s;  
  entry->set_hash();
  STR_HASH_ENTRY<HMM*> *res = tbl->get(entry);
  if(res != NULL)
    return res->val;
  else
    return (HMM*)NULL;
}
#endif
