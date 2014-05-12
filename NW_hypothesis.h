// -*- Mode: C++;  -*-
// File: NW_hypothesis.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* REQUIRES: 
//*
//* HISTORY:
//*  Jan 27 12:55 1995 (sjr): Added support for state path decoding
//*  Apr 13 16:22 1994 (sjr): Changed to integer repn. of logprobs
//* Created: Tues Apr 5, 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// First we need some LM-related global vars
GLOBAL_VAR_INIT(BOOL, State_decode, FALSE);
GLOBAL_VAR_INIT(BOOL, Pause_decode, FALSE);
GLOBAL_VAR_INIT(INT, Silence_index, -1);
GLOBAL_VAR_INIT(INT, Sentence_start_index, -1);
GLOBAL_VAR_INIT(INT, Sentence_end_index, -1);
GLOBAL_VAR_INIT(INT, Paragraph_index, -1);
GLOBAL_VAR_INIT(INT, Article_index, -1);
GLOBAL_VAR_INIT(INT, Unknown_word_index, -1);
GLOBAL_VAR(STRING, Silence_symbol);
GLOBAL_VAR(STRING, Sentence_start_symbol);
GLOBAL_VAR(STRING, Sentence_end_symbol);
GLOBAL_VAR(STRING, Paragraph_symbol);
GLOBAL_VAR(STRING, Article_symbol);
GLOBAL_VAR(STRING, Unknown_word_symbol);
GLOBAL_VAR(BOOL, Cross_sentence);

class HYP {
  static HYP* free_list;
  static INT block_size;
  static INT total_hyps;
  static LIST<INT>* gamma;

  
  // this does double service
  //  1)  in a constructed HYP it points to the previous element in the 
  //       hypothesis
  //  2)  in the memory pool it points to the next element in the free list
  HYP* previous;

  // this is used for online HYP output
  LIST<HYP*>* next;

  //  INT n_next;

  INT t;
  INT logprob;
  //  the difference between this logprob and the best of any hypothesis
  //   ending at this time
  INT deltalogprob;
  INT wd;
  INT pron_version;
  INT pron_prob;
  INT wd_logprob;
  INT* lm_context;
  INT count;

  // for State_decode
  QUEUE<INT>* nodeseq;
  QUEUE<INT>* nodetimeseq;
  QUEUE<INT>* nodeprobseq;

  //for mixture LMs
  LIST<INT>* ng_wts;

  //  LIST<INT>* next_lm_wts;
  //  LIST<INT>* current_lm_wts;

public:

  // redefine new and delete so that memory objects are allocated from and 
  //  deallocated to an explicitly defined free list
  void* operator new(size_t size);
  void operator delete(void* dead, size_t size);

  HYP(INT w, INT tim, INT totalp, INT wp, HYP *h, INT oov2unk_of_w = -1, 
      INT v = 0, INT r = 0){    
    total_hyps++;
    next = new LIST<HYP*>();
    previous = h;
    if(previous != NULL)
      previous->inc_next(this);

    //    n_next = 0;
    wd = w;
    t = tim;
    logprob = totalp;
    wd_logprob = wp;
    deltalogprob = 0;
    pron_version = v;
    pron_prob = r;

    if(oov2unk_of_w == -1)
      oov2unk_of_w = w;

    //  Now done by HYP::new
    //    lm_context = new INT[lm_context_length];
    if(previous == NULL) {
      if(w == Silence_index)
	lm_context[0] = Sentence_start_index;
      else
	lm_context[0] = oov2unk_of_w;

      // set default context
      for(INT i = 1; i < lm_context_length; i++)
	lm_context[i] = Sentence_start_index;
      if(lm_context[0] !=  Sentence_start_index)
	count = 1;
      else
	count = 0;

    } else{
      // copy and update previous context
      if(w == Silence_index) {
	// no change to LM context
	memcpy(lm_context, previous->get_lm_context(), 
	       lm_context_length*sizeof(INT));
      } else {
	memcpy(lm_context+1, previous->get_lm_context(), 
	       (lm_context_length-1)*sizeof(INT));
	lm_context[0] = oov2unk_of_w;
	count = previous->get_count() + 1;
      }
    }


    //    if(Debug_flag.mixlm)
    //      printf("count = %d, mixture_lm_size = %d\n", count, mixture_lm_size);

    if(mixture_lm_size > 0) {
      assert(ng_wts != NULL);

      ng_wts->set(logcode(1.0/(FLOAT)mixture_lm_size), mixture_lm_size);	
      if(count > 1)
    	update_ng_wts();
    }
  }

  ~HYP() { 
    total_hyps--;

    while(next->size > 0)
      next->pop()->unset_previous();
    delete next;

    if(previous != NULL)
      previous->dec_next(this); 

    //  Recycled
    //    delete lm_context;

    if(State_decode) {
      // Recycled
      //      delete nodeseq;
      //      delete nodetimeseq;
      //      delete nodeprobseq;
      nodeseq->fast_clear();
      nodetimeseq->fast_clear();
      nodeprobseq->fast_clear();
    }
  }
    
  static INT lm_context_length;
  static INT mixture_lm_size;
  static LIST<INT>* ng_probs;
  static INT lm_prob;
  static BOOL output_pron_version;
  static BOOL output_silence;

  INT get_total_hyps() { return total_hyps; }
  INT get_word() { return wd; }
  INT get_pron_version() { return pron_version; }
  INT get_pron_prob() { return pron_prob; }
  INT* get_lm_context() { return lm_context; }
  VOID set_lm_context(INT i, INT x) { lm_context[i] = x; }
  INT get_lm_context_length() { return lm_context_length; }
  VOID set_lm_context_length(INT n) { lm_context_length = n; }
  INT get_prev_word() { return get_lm_context()[1]; }
  LIST<HYP*>* get_next() { return next; }
  INT get_n_next() { return next->size; }
  VOID inc_next(HYP *h) { next->push(h); }
  VOID dec_next(HYP *h) { next->index_remove(h); }
  HYP* get_previous() { return previous; }
  VOID unset_previous() { previous = NULL; }
  INT get_logprob() { return logprob; }
  INT get_word_logprob() { return wd_logprob; }
  INT get_time() { return(t); }
  INT get_delta() { return(deltalogprob); }
  INT get_count() { return(count); }
  VOID set_delta(INT p) { deltalogprob = p; }
  VOID add_to_delta(INT p) { deltalogprob += p; }
  VOID compute_delta(INT p) { deltalogprob = p - logprob; }
  VOID set_logprob(INT s) { logprob = s; }
  VOID add_to_logprob(INT p) { logprob = logadd(logprob, p); }
  BOOL is_less_than(HYP* h) { return(logprob < h->logprob); }
  BOOL is_greater_than(HYP* h) { return(logprob > h->logprob); }
  BOOL is_equal(HYP* h) { return(logprob == h->logprob); }

  QUEUE<INT>* get_nodeseq() { return(nodeseq); }
  QUEUE<INT>* get_nodetimeseq() { return(nodetimeseq); }
  QUEUE<INT>* get_nodeprobseq() { return(nodeprobseq); }

  // mixture LM
  LIST<INT>* get_ng_wts() { return ng_wts; }
  VOID update_ng_wts();
  INT wt_sum_probs();



  VOID insert_pause(INT tim, INT p){
    INT wp = p - previous->get_logprob();
    HYP* h = new HYP(Silence_index, tim, p, wp, previous);
    previous->dec_next(this);
    previous = h;
    h->inc_next(this);
    wd_logprob -= wp;
  }

  // display methods
  VOID display(LIST<STRING> *vocab, GENFILE *fp, INT sc = -1,
	       BOOL op_logprob = FALSE, STRING first_char = "");
  VOID log_display(LIST<STRING> *vocab=NULL, GENFILE *fp=gstdout);
  VOID state_display(LIST<STRING> *phnms, INT offset, FLOAT fshift, 
		     GENFILE *fp);
  VOID ctm_display(LIST<STRING> *vocab, INT offset, FLOAT fshift,
		   GENFILE *fp, INT sc = -1, CHAR channel = 'A');
  VOID srt_display(LIST<STRING> *vocab, INT offset, FLOAT fshift, 
		   GENFILE *fp);

  VOID online_display(LIST<STRING> *vocab, GENFILE *fp){ 
    INT w = get_word();
    if(w == Sentence_start_index || w == Sentence_end_index || w == Silence_index)
      return;
    fp->fprintf("%s ", vocab->get(w&WORD_MASK)); 
    fp->fflush();
  }

  VOID online_log_display(LIST<STRING> *vocab=NULL, GENFILE *fp=gstdout) {
    INT w = get_word()&WORD_MASK;
    if(vocab == NULL) {
      if(output_pron_version)
	fp->fprintf("%d(%d) <%d> ", w, get_pron_version(), get_time());
      else
	fp->fprintf("%d <%d> ", w, get_time());
    } else {
      if(output_pron_version)
	fp->fprintf("%s/%d(%d) [%.2f] <%d> ", vocab->get(w), w,
		    get_pron_version(),
		    logdecodelog(get_logprob()), get_time());
      else
	fp->fprintf("%s/%d [%.2f] <%d> ", vocab->get(w), w,
		    logdecodelog(get_logprob()), get_time());
    }
    fp->fflush();
  }

  VOID online_state_display(LIST<STRING> *phnms, INT start, INT offset, 
			    FLOAT fshift, GENFILE *fp){
    INT w, end;
    LIST<INT>* statel = new LIST<INT>();
    LIST<INT>* timel = new LIST<INT>();    
    nodeseq->init_cursor();
    while(nodeseq->next())
      statel->push(nodeseq->read());
    nodetimeseq->init_cursor();
    while(nodetimeseq->next())
      timel->push(nodetimeseq->read());
    while(statel->size > 0) {
      w = statel->pop();
      end = timel->pop();
      if(start == end)
	continue;
      fp->fprintf("%.3f %.3f  %s (%d)\n", (start+offset) * fshift, 
		  (end-start) * fshift, phnms->get(w), w);
      start = end;
    }
    delete statel;
    delete timel;
    fp->fflush();
  }

  VOID online_ctm_display(LIST<STRING> *vocab, INT start, INT offset, FLOAT fshift,
			  GENFILE *fp, INT sc = -1, CHAR channel = 'A'){ 
    if(start == get_time()) return;
    INT w = get_word();
    if(w == Sentence_start_index || w == Sentence_end_index)
      w = Silence_index;
    if(w == Silence_index && !output_silence)
      return;
    w &= WORD_MASK;
    if(output_pron_version)
      fp->fprintf("%d %c %.3f %.3f %s(%d) %.4f\n", 
		  sc, channel, (start+offset)*fshift, 
		  (get_time()-start)*fshift, vocab->get(w), 
		  get_pron_version(),
		  logdecodelog(get_word_logprob())); 
    else
      fp->fprintf("%d %c %.3f %.3f %s %.4f\n", sc, channel, (start+offset)*fshift, 
		  (get_time()-start)*fshift, vocab->get(w), 
		  logdecodelog(get_word_logprob()));      
    fp->fflush();
  }

  VOID online_srt_display(LIST<STRING> *vocab, INT start, INT offset, 
			  FLOAT fshift, GENFILE *fp){
    if(start == get_time()) return;
    INT w = get_word();
    if(w == Sentence_start_index || w == Sentence_end_index)
      w = Silence_index;
    if(w == Silence_index && !output_silence)
      return;
    w &= WORD_MASK;
    if(output_pron_version)
      fp->fprintf("<Word S_time=%.3f E_time=%.3f Prob=%.4g Version=%d> %s </Word>\n", 
		  (start+offset)*fshift, (get_time()+offset)*fshift, 
		  logdecodelog(get_word_logprob()),
		  get_pron_version(),
		  vocab->get(w));
    else
      fp->fprintf("<Word S_time=%.3f E_time=%.3f Prob=%.4g> %s </Word>\n", 
		  (start+offset)*fshift, (get_time()+offset)*fshift, 
		  logdecodelog(get_word_logprob()),
		  vocab->get(w));      
    fp->fflush();
  }

  VOID info(LIST<STRING> *vocab = NULL, LIST<STRING> *phns = NULL, 
	    GENFILE *fp = gstdout, STRING hdr = "");
};

GLOBAL_VAR(LIST<HYP*>*, Hyp_garbage_list);

