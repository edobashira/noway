// -*- Mode: C++;  -*-
// File: hmm.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION: Basic HMM classes
//*
//* CLASSES: HMM_TRANS, HMM
//* 
//* REQUIRES: *INTERFACE classes (acoustic.h)
//*           LEXICON (lexicon,h) - friend class
//*
//* HISTORY:
//*  Feb 21 12:50 1995 (sjr): *_in() now include lub_n parm for init_lub()
//*  Jan 27 12:57 1995 (sjr): Added global state indexing and support for
//*                             state path decoding
//*  Apr 21 12:02 1994 (sjr): Tested okay using ASCII_FLOAT_INTERFACE.
//*  Apr 20 17:30 1994 (sjr): added shared and interface routines.
//*                           tested okay.
//*                           ACOUSTIC_INTERFACE stuff not tested.
//*  Apr  8 17:15 1994 (sjr): data structures only, tested okay
//* Created: Wed Apr  6 14:31:29 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


GLOBAL_VAR_INIT(BOOL,Backtrace, FALSE);
GLOBAL_VAR_INIT(BOOL,Forward_process, FALSE);
GLOBAL_VAR_INIT(INT, Silence_distr, 0);

class TREE_LEXICON;

class HMM_TRANS {
public:
  HMM_TRANS(INT f = 0, INT t = 0, INT p = 0) {
    from = f;
    to = t;
    prob = p;
  }
  
  ~HMM_TRANS(){}

  INT from;
  INT to;
  INT prob;

  BOOL equal(HMM_TRANS* t) {
    return(from == t->from && to == t->to && prob > INT_MIN);
  }

  BOOL equal(INT f, INT t) {
    return(from == f && to == t && prob > INT_MIN);
  }


};

const INT DEF_BACKTRACE_LIST_ALLOC = 32;
class BACKTRACE_LIST : public LIST<INT> {
  static BACKTRACE_LIST* free_list;
  static INT block_size;
  static INT total_backtrace_lists;
  
  BACKTRACE_LIST* previous;
public:
  void* operator new(size_t size);
  void operator delete(void* dead, size_t size);

  BACKTRACE_LIST() : LIST<INT>(-1) {
    previous = NULL;
  }

  ~BACKTRACE_LIST() {}

  VOID free();

  INT get_total_backtrace_lists() { return total_backtrace_lists; }
};




class HMM : public ACOUSTIC {
  friend class ACOUSTIC_LEXICON;
  friend class WORD;
protected:
  static ACOUSTIC_INTERFACE *input;

  static INT nmodels;
  static INT ndistr;
  static INT total_states;	// doesn't include start and finish states

  INT nstates;			// includes start and finish states
  INT id;
  INT distr;
  INT *logprob;
  INT *newlogprob;
  INT *statebp;
  INT *newstatebp;
  INT ntrans;
  //  INT *trans_from;
  //  INT *trans_to;
  //  INT *trans_prob;
  HMM_TRANS *trans;

  // used in state decoding
  // stuff for strict lub
  INT maxprob;
  BACKTRACE_LIST* parent_exit_times;
  BACKTRACE_LIST* maxprobs;
  BACKTRACE_LIST* exitprobs;
  INT entry_time;
  INT frame_activated;
  INT last_activated;

  FLOAT exit_prob;
  FLOAT self_prob;

public:
  HMM(INT ns = 5, STRING nm = "", INT di = -1, BOOL make_trans = TRUE) {
    if(nm != NULL) 
      name = strdup(nm);
    else
      name = strdup("");
    id = nmodels;
    if(di == -1)
      distr = id;
    else
      distr = di;
    nmodels++;
    total_states = total_states + (ns - FIRST_STATE);
    nstates = ns;
    logprob = new INT[nstates];
    statebp = new INT[nstates];
    clear_logprob();
    last_activated = frame_activated = -1;
    newlogprob = new INT[nstates];
    newstatebp = new INT[nstates];
      memset(newlogprob, 0x7f, nstates*sizeof(INT));
      memset(newstatebp, 0x7f, nstates*sizeof(INT));
    if(make_trans) {
      ntrans = 0;
      //      trans_from = new INT[nstates*nstates];
      //      trans_to = new INT[nstates*nstates];
      //      trans_prob = new INT[nstates*nstates];
      trans = new HMM_TRANS[nstates*nstates];
    }
//    if(Backtrace) {
//      parent_exit_times = new LIST<INT>(8);
//      maxprobs= new LIST<INT>(8);
//      if(State_decode || Pause_decode)
//	exitprobs = new LIST<INT>(8);
//    }
    parent_exit_times = NULL;
    maxprobs = NULL;
    exitprobs = NULL; 

    self_prob = 0.5;
    exit_prob = 0.5;
  }

  ~HMM() {}

  VOID phi_init(FLOAT self_prob, FLOAT exit_prob, FLOAT dur_scale, FLOAT pdp);

  STRING name;
  INT get_nmodels() { return(nmodels); }
  INT get_nstates() { return(nstates); }
  INT get_id() { return(id); }
  INT get_sentence_length() { return(sentence_length); }
  INT get_distr() { return(distr); }
  BOOL reached_exit() { return(logprob[EXIT_STATE] != PROB_NULL); }
  INT get_entry_prob() { return(logprob[ENTRY_STATE]); }
  INT get_exit_prob() { return(logprob[EXIT_STATE]); }
  INT get_first_prob() { return(logprob[FIRST_STATE]); }
  INT get_state_prob(INT i) { return(logprob[i]); }
  INT* get_newlogprob() { return(newlogprob); }
  INT* get_newstatebp() { return(newstatebp); }
  //  INT* get_trans_prob() { return(transprob); }
  //  INT* get_trans_to() { return(trans_to); }
  //  INT* get_trans_from() { return(trans_from); }
  HMM_TRANS* get_trans() { return(trans); }
  INT get_ntrans() { return(ntrans); }

  VOID set_init_logprob(INT p, INT t) { 
    logprob[ENTRY_STATE] = p; 
    statebp[ENTRY_STATE] = t;
  }

  VOID set_startprob(INT p) { startprob = p; }
  VOID add_init_logprob(INT p) { 
    logprob[ENTRY_STATE] = logadd(logprob[ENTRY_STATE],p); 
  }
  VOID clear_logprob() { 
    memset(logprob, 0x7f, nstates*sizeof(INT)); 
    if(Backtrace)
      memset(statebp, 0x7f, nstates*sizeof(INT)); 
  }

  VOID hmm_set_lub(INT t);

  INT compute_maxprob() {
    //    if(logprob[ENTRY_STATE] == PROB_NULL)
    maxprob = INT_MIN;
    //else
    //maxprob = logprob[ENTRY_STATE];

    for(INT i = FIRST_STATE;  i < nstates; i++) {
      if(logprob[i] != PROB_NULL)
	maxprob = MAX(logprob[i], maxprob);
    }
    return(maxprob);
  }

  VOID set_maxprob(INT p) { maxprob = p; }
  VOID inc_maxprob(INT p) { maxprob += p; }
  INT get_maxprob() { return(maxprob); }

  VOID add_trans(INT from_state, INT to_state, FLOAT prob) {
    INT lp = INT_MIN;
    if(prob > 0.0)
      lp = logcode(prob);
    if(lp < logtrans && lp > INT_MIN) logtrans = lp;
    trans[ntrans].from = from_state;
    trans[ntrans].to = to_state;
    trans[ntrans].prob = lp;
    ntrans++;
  }

  BOOL check_trans(HMM_TRANS* tr) {
    for(INT i = 0; i < ntrans; i++)
      if(trans[i].equal(tr)) return(TRUE);
    return(FALSE);
  }


  BOOL check_trans(INT fs, INT ts) {
    for(INT i = 0; i < ntrans; i++)
      if(trans[i].equal(fs, ts)) return(TRUE);
    return(FALSE);
  }

  BOOL step_forward(INT t);

  // return the difference between prob p and the least upper bound so far
  // at time t
  INT delta(INT p, INT t) {
    return(lub->get(t) - p);
  }

  VOID update_lub(INT p, INT t) { 
    if(Debug_flag.lub)
      printf("%s/%d (%d): lub[%d] = %d (was %d)\n", 
	     name, id, start_time, t, p, lub->get(t));
    lub->put(t, p);
    if(t+1 < proc_frames_seen()) {
      INT tt = t+1;
      INT pp = p + def_lub_inc->get(t+1);
      if(lub->get(tt) < pp){
	if(Debug_flag.lub)
	  printf(" %s/%d (%d): lub[%d] = %d (was %d)\n", 
		 name, id, start_time, tt, pp, lub->get(tt));
	lub->put(tt,pp);
      }
    }
  }

  // get the output prob for this model at time t, returning 0 if
  //  t is past end of sentence
  INT output_prob(INT t) { 
    if(t < proc_frames_seen())
      return(outputs->get(t)[distr]); // output prob in cache
    else if(t >= sentence_length)
      return(0);		// past end of sentence
    else {
      INT* localp = new INT[input->get_nout()];
      if(input->read(localp, t) == FALSE){ // reached end of sentence
	delete[] localp;
	sentence_length = t;
	return(0);
      } else {  // read probs from input stream
	outputs->enqueue(localp);
	inc_frames_seen();
	inc_proc_frames_seen();
	return(localp[distr]);
      }
    }
  }

  VOID skip_sentence_acoustics() {
    new_sentence();
    skip_frames_in_current_sent(INT_MAX);
    sentence_count++;
    return;
  }

  INT skip_frames_in_current_sent(INT n) {
    INT t = 0;
    while(t < n && (ndx_output || input->read(NULL, t) != FALSE)){
      inc_frames_seen();
      t++;
    }
    return(n-t);
  }
    

  BOOL is_floor(INT p) {
    return(p == input->get_floor());
  }

  VOID ascii_float_in(LIST<FLOAT> *priors, STRING floor_file, 
		      STRING probs_file, 
		      FLOAT sc, BOOL lf, INT nfl) {
    amscale = sc;
    if(input != NULL)
      fprintf(stderr, "WARNING: HMM::ascii_float_in resetting acoustic input\n");
    input = new ASCII_FLOAT_INTERFACE();
    input->initialize(priors, floor_file, probs_file, sc, lf, nfl);
  }

#ifdef ASCII_HEX
  VOID ascii_hex_in(LIST<FLOAT> *priors, STRING floor_file,
		    STRING probs_file, 
		    FLOAT sc, BOOL lf, INT nfl){
    amscale = sc;
    if(input != NULL)
      fprintf(stderr, "WARNING: HMM::ascii_hex_in resetting acoustic input\n");
    input = new ASCII_HEX_INTERFACE;
    input->initialize(priors, floor_file, probs_file, sc, lf, nfl);
  }
#endif

  VOID bin_float_in(LIST<FLOAT> *priors, STRING floor_file,
		    STRING probs_file, 
		    FLOAT sc, BOOL lf, INT nfl){
    amscale = sc;
    if(input != NULL)
      fprintf(stderr, "WARNING: HMM::bin_float_in resetting acoustic input\n");
    input = new BIN_FLOAT_INTERFACE;
    input->initialize(priors, floor_file, probs_file, sc, lf, nfl);
  }


  VOID lnafloat_in(LIST<FLOAT> *priors, STRING floor_file,
		    STRING probs_file, 
		    FLOAT sc, BOOL lf, INT nfl, INT no){
    amscale = sc;
    if(input != NULL)
      fprintf(stderr, "WARNING: HMM::lnafloat_in resetting acoustic input\n");

    input = new LNAFLOAT_INTERFACE;
    input->initialize(priors, floor_file, probs_file, sc, lf, nfl, no);
  }


  VOID lna_in(LIST<FLOAT> *priors, STRING floor_file, STRING probs_file, 
	      FLOAT sc, BOOL lf, INT nfl, INT no){
    amscale = sc;
    if(input != NULL)
      fprintf(stderr, "WARNING: HMM::lna_in resetting acoustic input\n");
    input = new LNA_INTERFACE;
    input->initialize(priors, floor_file, probs_file, sc, lf, nfl, no);
  }

#ifdef SOCKETIO
  VOID socket_in(LIST<FLOAT> *priors, STRING floor_file, INT port,
		 FLOAT sc, BOOL lf, INT nfl){
    amscale = sc;
    if(input != NULL)
      fprintf(stderr, "WARNING: HMM::socket_in resetting acoustic input\n");
    input = new SOCKET_INTERFACE;
    char buf[32];
    sprintf(buf, "%d", port);
    input->initialize(priors, floor_file, buf, sc, lf, nfl);
  }
#endif

  VOID set_sum_check() { input->set_sum_check(); }
  VOID unset_sum_check() { input->unset_sum_check(); }

  VOID backtrace_garbage_collect() {
    //    delete parent_exit_times;
    //    delete maxprobs;
    //    if(Pause_decode || State_decode) 
    //      delete exitprobs;
    parent_exit_times->free();
    maxprobs->free();
    if(Pause_decode || State_decode) 
      exitprobs->free();
  }

  VOID print() {
    printf("name: %s\n", name);
    printf("nstates: %d\n", nstates);
    printf("id: %d\n", id);
    printf("distr: %d\n", distr);
    printf("nmodels: %d\n", nmodels);
    //    printf("n_lub: %d\n", lub->size);
    //    printf("sentence_length: %d\n", sentence_length);
    //    printf("input: 0x%p\n", input);
//    print_trans();
    //    printf("n_outputs: %d\n", outputs->size);
//    print_outputs();
    printf("\n");
  }

  VOID print_trans() {
  }

  VOID print_outputs(){
    INT i = 0;
    while(i < outputs->size){
      printf("%.2f ", logdecode(outputs->get(i)[distr]));
      if(!(++i%10)) printf("\n");
    }
    if(outputs->size > 0) printf("\n");
  }

};


