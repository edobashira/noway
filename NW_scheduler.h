// -*- Mode: C++;  -*-
// File: scheduler.h
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
//* Created: Fri Apr 22 10:58:32 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef TIMING
EXTERN INT tick;
EXTERN INT sentence_tick;
EXTERN struct tms start_time_buffer;
EXTERN struct tms end_time_buffer;
EXTERN clock_t startt;
EXTERN clock_t endt;
EXTERN time_t tt;
#endif

GLOBAL_VAR_INIT(INT, Version, 2);
GLOBAL_VAR_INIT(INT, Subversion, 9);
GLOBAL_VAR_INIT(INT, Patchlevel, 4);

class SECTION {
  INT stime;
  FLOAT fstime;
  INT etime;
  FLOAT fetime;
  STRING type;
  STRING topic;
  STRING id;
public:
  SECTION(FLOAT s, FLOAT e, FLOAT shift, STRING ty, STRING to, STRING i) {
    FLOAT x = shift*0.1;
    fstime = s;
    fetime = e;
    stime = (INT)((s+x)/shift);
    etime = (INT)((e+x)/shift);
    if(Debug_flag.ndx) 
      printf("e = %f   x = %f   shift = %f  e+x/shift = %f INT = %d\n",
	     e, x, shift, (e+x)/shift, etime);
    type = strdup(ty);
    topic = strdup(to);
    id = strdup(i);
  }

  INT get_start_time() { return stime; }
  INT get_end_time() { return etime; }
  FLOAT get_float_start_time() { return fstime; }
  FLOAT get_float_end_time() { return fetime; }
  STRING get_type() { return type; }
  STRING get_topic() { return topic; }
  STRING get_id() { return id; }

  VOID print_start(GENFILE* fp, CHAR* comment = NULL){ 
    if(comment == NULL)
      fp->fprintf("<Section S_time=%.3f E_time=%.3f Type=%s Topic=%s ID=%s >\n",
		  fstime, fetime, type, topic, id); 
    else
      fp->fprintf("%s <Section S_time=%.3f E_time=%.3f Type=%s Topic=%s ID=%s >\n",
		  comment, fstime, fetime, type, topic, id); 
  }

  VOID print_end(GENFILE* fp, CHAR* comment = NULL){ 
    if(comment == NULL)
      fp->fprintf("</Section>\n"); 
    else
      fp->fprintf("%s </Section>\n", comment);       
  }

  VOID info(GENFILE* fp = gstdout) {
    fp->fprintf("ID = %s\n", id);
    fp->fprintf("Start_time = %d (%f)\n", stime, fstime);
    fp->fprintf("End_time = %d (%f)\n", etime, fetime);
    fp->fprintf("Type = %s\nTopic = %s\n\n", type, topic);
  }
};

class SCHEDULER : public ACOUSTIC {
 private:
  // A separate priority queue of hypotheses at each time
  QUEUE <STACK*> *stacks;

  // pronunciation dictionary and language model
  ACOUSTIC_LEXICON* dictionary;
  GRAMMAR* language;

  PARAM_TABLE *param;
  static STRING decoder_param[];

  GENFILE* output_fp;
  GENFILE* state_dec_fp;
  BOOL output_logprob;
  STRING lattice_name;
  //  STRING detailed_lattice_name; 
  STRING ctm_filename;
  GENFILE* ctm_fp;
  STRING srt_filename;
  GENFILE* srt_fp;
  STRING ndx_hdr;
  STRING ndx_ftr;
  SECTION* ndx_sec;
  QUEUE<SECTION*>* ndxq;
  BOOL dump_lub;
  STRING dump_lub_name;
  INT inc_output;
  INT maxhyps;
  INT n_hyps_db;

  BOOL online;
  BOOL online_output;
  HYP *firsth;
  INT firsth_start;

  QUEUE<INT> *depth_count;
  ULONG total_nodes_used;
  ULONG total_nodes_used2;
  ULONG total_hyps_created;
  ULONG total_hyps_used;
  ULONG total_phones_pruned;
  ULONG total_nodes_pw;
  INT max_nhyps_created; 
  INT max_nnodes_used;
  INT max_nnodes_used2;

  LIST<LIST<LATTICE_LINK*>*> *lattices;
  LIST<LIST<INT>*> *silends;
  BOOL build_lattice;
  BOOL pipe_lattice; // DCA 18/MAR/96
  INT nlatt_nodes;

  BOOL state_decode_output;
  BOOL ctm_output; // DCA 25/MAR/96
  BOOL srt_output;
  BOOL inc_offset;
  INT frame_offset;
  INT sentence_frame_offset;
  FLOAT frame_shift;

  BOOL cache_lm;

  BOOL maintain_article_context;

 public:
  SCHEDULER() {
    dictionary = NULL;
    language = NULL;
    param = NULL;
    stacks = NULL;
    output_fp = gstdout;
    n_hyps_db = 5;
    inc_output = 0;
    firsth = NULL;
    firsth_start = 0;
    online = FALSE;
    online_output = FALSE;
    output_logprob = FALSE;
    depth_count  = new QUEUE<INT>();
    total_hyps_created = 0;
    total_hyps_used = 0;
    total_nodes_used = 0;
    total_nodes_used2 = 0;
    total_phones_pruned = 0;
    total_nodes_pw = 0;
    max_nhyps_created = 0;
    max_nnodes_used = 0;
    max_nnodes_used2 = 0;
    lattices = new  LIST<LIST<LATTICE_LINK*>*>(2048);
    silends = new  LIST<LIST<INT>*>(2048);
    build_lattice = FALSE;
    pipe_lattice = FALSE; 
    state_decode_output = FALSE;
    ctm_output = FALSE; // DCA 25/MAR/96
    srt_output = FALSE;
    inc_offset = FALSE;
    ctm_filename = NULL;
    srt_filename = NULL;
    ndx_hdr = NULL;
    ndx_ftr = NULL;
    ndx_sec = NULL;
    ndxq = NULL;
    ctm_fp = NULL;
    srt_fp = NULL;
    sentence_frame_offset = 0;
    frame_offset = 0;
    frame_shift = 0.016; 
    nlatt_nodes = 0;
    cache_lm = FALSE;
  }

  ~SCHEDULER() {
  }

  VOID main_prog(INT argc, STRING argv[]);

  // returns a priority queue of decoding hypotheses for the current sentence
  STACK* decode();

  VOID welcome_mesg() { 
    if(Patchlevel == 0)
      printf("noway version %d.%d\n", Version, Subversion); 
    else
      printf("noway version %d.%dp%d\n", Version, Subversion, Patchlevel); 
    fflush(stdout); 
  }

  VOID parse(INT argc, STRING argv[]);
  VOID build();

  VOID garbage_collect(LIST<HYP*>* hlist);
  VOID hyp_chain_forward();
  VOID init_ndx(GENFILE* ndx_fp);

  // Print out best hypothesis
  VOID init_online_output();
  VOID output_hyps(STACK *hyps);

  // DCA 14/MAR/96: Print out lattice header
  VOID print_lattice_header(GENFILE *fp, PARAM_TABLE *param, HYP* hyp, 
			    LIST<STRING> *vocab, INT sentence_count, 
			    INT num_links, INT num_nodes);

};

class WORD_TIME_PAIR {
  INT wd;
  INT t;
public:
  WORD_TIME_PAIR(INT word, INT tim) {
    wd = word;
    t = tim;
  }

  ~WORD_TIME_PAIR() {}

  INT get_word() { return wd; }
  INT get_time() { return t; }
  VOID set_time(INT tim) { t = tim; }
  VOID set_word(INT word) { wd = word; }

  friend BOOL operator==(WORD_TIME_PAIR &x, WORD_TIME_PAIR &y) {
    return((x.wd == y.wd) && (x.t == y.t)); 
  }
  friend BOOL operator!=(WORD_TIME_PAIR &x,  WORD_TIME_PAIR &y) {
    return((x.wd != y.wd) || (x.t != y.t)); 
  }

};
