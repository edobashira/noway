// -*- Mode: C++;  -*-

class SUBWORD : public HMM {
protected:
  static LIST<SUBWORD*>* active_sw_list;
  static INT nsubwords;
public:

  SUBWORD(HMM *p) : HMM(p->get_nstates(), p->name, p->get_distr(), FALSE) {  
      trans = p->get_trans();
      ntrans = p->get_ntrans();
      deactivate();
      logprob[EXIT_STATE] = PROB_NULL;
      nsubwords++;
  }

  LIST<SUBWORD*>* get_active_sw_list() 
    { return active_sw_list; }

  BOOL activated;

  virtual VOID node_sequence(INT, INT, HYP*) {}
  virtual VOID backtrace_lub(INT, INT) {}

  BOOL forward(INT t);
  BOOL prune(INT t);

  virtual BOOL activate(INT t, INT p) { 
    BOOL react = FALSE;
    set_init_logprob(p, t);
    activated = TRUE;
    if(Backtrace) 
      memset(&statebp[1], 0x7f, (nstates-1)*sizeof(INT));

    // nodes can be reactivated with the same start time
    if(frame_activated == start_time && last_activated == sentence_count) {
      react = TRUE;
      if(Backtrace) {
	INT gap = (t - entry_time) - maxprobs->size;
	if(gap > 0) {
	  maxprobs->set(PROB_NULL, gap, maxprobs->size);
	  parent_exit_times->set(-1, gap, parent_exit_times->size);
	  if(State_decode || Pause_decode){
	    exitprobs->set(PROB_NULL, gap, exitprobs->size);
	  }
	}
      }
    } else {
      entry_time = t;
      frame_activated = start_time;
      last_activated = sentence_count;
      if(Backtrace) {
	parent_exit_times = new BACKTRACE_LIST();
	maxprobs = new BACKTRACE_LIST();
	if(Pause_decode || State_decode) 
	  exitprobs = new BACKTRACE_LIST();
	active_sw_list->push(this);
      }
    }
    return(react);
  }
  
  VOID deactivate() {
    activated = FALSE;
    // Reset all states except exit state
    logprob[ENTRY_STATE] = PROB_NULL;
    memset(&logprob[FIRST_STATE], 0x7f, (nstates-2)*sizeof(INT));    
  }

  virtual VOID short_info(INT) {}
  virtual VOID info(INT, STRING) {}
  
};


