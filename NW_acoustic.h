// -*- Mode: C++;  -*-
// File: acoustic.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Classes for reading in acoustic information,
//*             typically phone probabilities
//*
//* CLASSES:  ACOUSTIC, ACOUSTIC_INTERFACE, FILE_INTERFACE,
//*           ASCII_FLOAT_INTERFACE, ASCII_HEX_INTERFACE,
//*           BIN_FLOAT_INTERFACE, LNA_INTERFACE, SOCKET_INTERFACE
//* 
//* REQUIRES: 
//*
//* HISTORY:
//*  Jun 19 17:20 1995 (sjr): #ifdefed ASCII_HEX_INTERFACE
//*  Feb 21 12:49 1995 (sjr): Extra parm in initialise (nfl) for init_lub
//*  Jul 14 12:02 1994 (sjr): ACOUSTIC class used for statics;  also
//*                            inherited by HMM and SCHEDULER
//*  Apr 21 12:01 1994 (sjr): Tesed ASCII_FLOAT_INTERFACE okay.
//*  Apr 19 15:31 1994 (sjr): Implemented interfaces from y0 
//*                           (ascii float/hex, bin float, socket)
//*                           Not yet done lna.  Untested.
//* Created: Tue Apr 19 11:18:55 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "NW_collection.h"

class ACOUSTIC {
 protected:
  static INT nout;
  static QUEUE_OFFSET<INT*> *outputs;
  static QUEUE_OFFSET<INT>  *lub;	// least upper bound used in pruning
  static QUEUE_OFFSET<INT>  *def_lub_inc;	
  static INT frame_count;
  static INT frames_this_sentence;
  static INT proc_frame_count;
  static INT proc_frames_this_sentence;
  static INT  sentence_length;
  VOID set_sentence_length(INT n) { sentence_length = n; }
  static INT sentence_count;
  static INT logtrans;

  // used in pruning
  static INT pathword_threshold;// below this num look at LM when pruning this node
  static INT beam;	  	// log likelihood beam for pruning
  static INT state_beam;	// log likelihood beam for pruning
  static INT nnodes_used;
  static INT nphones_pruned;
  static FLOAT *prob_accum;
  static FLOAT *data_priors;
  static INT start_time;

  static FLOAT prob_min;
  static FLOAT lna_floor;
  static INT logcode_min;

  static BOOL strict_lub;

  static INT startprob;

  static INT garbage_model; 

  //  static INT nreactivated;

  static INT average_lm_contexts;


  static INT lmscale;
  static INT lmoffset;
  static FLOAT amscale;
  static FLOAT logbase;

  static BOOL ndx_output;


public:

  // clean up the shared variables before starting a new sentence
  VOID new_sentence(BOOL reset_outputs = TRUE, BOOL reset_slen = TRUE) {
    nphones_pruned = 0;
    //    nreactivated = 0;
    if(reset_outputs) {
      while(outputs->size > 0)
	delete outputs->dequeue();
      outputs->fast_clear();
      lub->fast_clear();
      def_lub_inc->fast_clear();
      frames_this_sentence = 0;
      proc_frames_this_sentence = 0;
      if(reset_slen)
	set_sentence_length(INT_MAX);
    }
  }

  INT frames_seen() { return(frames_this_sentence); }
  INT total_frames_seen() { return(frame_count); }
  INT proc_frames_seen() { return(proc_frames_this_sentence); }
  INT total_proc_frames_seen() { return(proc_frame_count); }
  VOID inc_frames_seen() { frames_this_sentence++; frame_count++; }
  VOID inc_proc_frames_seen() { proc_frames_this_sentence++; proc_frame_count++; }


  INT get_beam() { return beam; }
  INT get_state_beam() { return state_beam; }
  FLOAT get_lna_floor() { return lna_floor; }

};

class ACOUSTIC_INTERFACE : public ACOUSTIC {
 protected:
  INT lognout;
  FLOAT *probs;
  FLOAT prob_sum;
  FLOAT sum_tolerance;
  BOOL sum_check;

  INT *logpriors;
  FLOAT scale;

  BOOL lower_floor;
  FLOAT *phone_threshold;
  FLOAT floor_prob;
  INT floor_logcode;

  VOID init_lub(INT frame);
  INT n_for_lub_init;
 public:
  ACOUSTIC_INTERFACE() { 
    scale = 1.0; 
    floor_prob = 1.0;
    floor_logcode = 0;
    lower_floor = FALSE;
    n_for_lub_init = 5;
    sum_check = TRUE; 
    sum_tolerance = 0.2;
  }

  VOID set_sum_check() { 
    sum_check = TRUE; 
    sum_tolerance = 0.2;
    if(Verbose)
      printf("Checking lna data sums to 1 +/- %g\n", sum_tolerance);
  }

  VOID unset_sum_check() { 
    if(Verbose)
      printf("NOT checking lna data sums to 1 +/- %g\n", sum_tolerance);
    sum_check = FALSE; 
  }

  INT get_nout() {return(nout); }
  INT get_floor() { return(floor_logcode); }
  virtual VOID initialize(LIST<FLOAT> *priors, STRING floor_name, 
			  STRING fname, 
			  FLOAT sc = 1.0, BOOL lf = FALSE, 
			  INT nfl = 5, INT no = 79) = 0;
  virtual BOOL read(INT *buffer, INT frame) = 0;
};

class FILE_INTERFACE : public ACOUSTIC_INTERFACE {
 protected:
  GENFILE *fp;
 public:
  FILE_INTERFACE() : ACOUSTIC_INTERFACE() { fp = NULL; }
  virtual VOID initialize(LIST<FLOAT> *priors, STRING floor_name = NULL, 
			  STRING fname = "-", 
			  FLOAT sc = 1.0, BOOL lf = FALSE, 
			  INT nfl = 5, INT no = 79);
  virtual BOOL read(INT *buffer, INT frame) = 0;
};

class ASCII_FLOAT_INTERFACE : public FILE_INTERFACE {
 public:
  ASCII_FLOAT_INTERFACE() : FILE_INTERFACE() {}
  BOOL read(INT *buffer, INT frame);
};

#ifdef ASCII_HEX
class ASCII_HEX_INTERFACE : public FILE_INTERFACE {
 public:
  ASCII_HEX_INTERFACE() : FILE_INTERFACE() {}
  BOOL read(INT *buffer, INT frame);
};
#endif

class BIN_FLOAT_INTERFACE : public FILE_INTERFACE {
 public:
  BIN_FLOAT_INTERFACE() : FILE_INTERFACE() {}
  VOID initialize(LIST<FLOAT> *priors, STRING floor_name = NULL, 
		  STRING fname = "-", 
		  FLOAT sc = 1.0, BOOL lf = FALSE, 
		  INT nfl = 5, INT no = 79);
  BOOL read(INT *buffer, INT frame);
};

class LNAFLOAT_INTERFACE : public FILE_INTERFACE {
 public:
  LNAFLOAT_INTERFACE() : FILE_INTERFACE() {}
  VOID initialize(LIST<FLOAT> *priors, STRING floor_name = NULL, 
		  STRING fname = "-", 
		  FLOAT sc = 1.0, BOOL lf = FALSE, 
		  INT nfl = 5, INT no = 79);
  BOOL read(INT *buffer, INT frame);
};

class LNA_INTERFACE : public FILE_INTERFACE {
protected:
  UCHAR* byte_buf;
 public:
  LNA_INTERFACE() : FILE_INTERFACE() {
    byte_buf = NULL;
  }
  VOID initialize(LIST<FLOAT> *priors, STRING floor_name = NULL, 
		  STRING fname = "-", 
		  FLOAT sc = 1.0, BOOL lf = FALSE, 
		  INT nfl = 5, INT no = 79);
  BOOL read(INT *buffer, INT frame);
};


#ifdef SOCKETIO
class SOCKET_INTERFACE : public ACOUSTIC_INTERFACE {
 protected:
  SOCKET *sock;
  UINT data_size;
 public:
  SOCKET_INTERFACE() : ACOUSTIC_INTERFACE() {}
  VOID initialize(LIST<FLOAT> *priors, STRING floor_name = NULL, 
		  STRING fname = "5765", 
		  FLOAT sc = 1.0, BOOL lf = FALSE, 
		  INT nfl = 5, INT no = 79);
  BOOL read(INT *buffer, INT frame);
};
#endif


