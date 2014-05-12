// -*- Mode: C++;  -*-
// File: acoustic.cc
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  Initialisation and reading routines for input
//*             probabilities
//*
//* CLASSES:   ACOUSTIC_INTERFACE and children
//* 
//* REQUIRES:
//*
//* HISTORY:
//*  Feb 21 12:49 1995 (sjr): Added support for variable number of probs to 
//*                            average over in init_lub()
//*  Jul 22 11:24 1994 (sjr): Implemented read + init for lna (single files).
//*                           Tested okay.
//*  Jul 14 12:00 1994 (sjr): Includes modifications for posterior-directed 
//*                            pruning
//*  Apr 21 12:01 1994 (sjr): Implemented read + init for float, hex, bin
//*                            and socket.  Tested float.
//* Created: Tue Apr 19 11:37:13 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "all.h"

static 
INT flt_comp(const VOID* f, const VOID* g) {
  FLOAT ff = *((FLOAT*)f);  
  FLOAT gg = *((FLOAT*)g);  
  if(ff==gg) 
    return(0);
  else
    return( ff < gg ? 1 : -1 );
}

VOID 
ACOUSTIC_INTERFACE::init_lub(INT frame) {

  // check for sum to 1 if required
  if(sum_check && fabs(prob_sum - 1.0) > sum_tolerance)
    panic("LNA file does not sum to 1 (%f on frame %d)\n", prob_sum, frame);

  if(proc_frames_seen() != frame)
    panic("proc_frames_seen() [%d] != frame [%d]", proc_frames_seen(), frame);
  

  // init lub to some garbage model estimate
  INT x;
  if(garbage_model >= 0) {
    x = logcode(probs[garbage_model]) - logpriors[garbage_model];
  } else {
    qsort(probs, nout, sizeof(FLOAT), flt_comp);
    FLOAT p = 0.0;
    for(INT i = 1; i < n_for_lub_init; i++){
      if(probs[i] != floor_prob)
	probs[i] = pow(probs[i], scale);
      p += probs[i];
    }
    x = logcode(p/(FLOAT)n_for_lub_init) - lognout;
  }
  x += logtrans;
  def_lub_inc->enqueue(x);
  if(frame > 0)
    x += lub->last();
  lub->enqueue(x);
  if(Debug_flag.lub)
    printf("lub[%d] = %d  *DEFAULT*\n", frame, x);
}

VOID
FILE_INTERFACE::initialize(LIST<FLOAT> *priors, STRING floor_name, 
			   STRING fname, FLOAT sc, BOOL lf, 
			   INT nfl, INT ) {
  CHAR buf[10000];
  GENFILE *pfp;
  
  scale = sc;
  lower_floor = lf;
  n_for_lub_init = nfl;
  if(scale != 1.0){
    floor_prob = pow(prob_min, scale);
    floor_logcode = (INT) ((FLOAT)logcode_min * scale + 0.5);
  } else {
    floor_prob = prob_min;
    floor_logcode = (INT)((FLOAT)logcode_min + 0.5);
  }
  if (Verbose) {
    printf("Initialising input probs(ascii)...\n");
  }

  if(!strcmp(fname, "stdin") || !strcmp(fname, ""))
    fp = gstdin;
  else
    fp = new GENFILE(fname, "r");

  do{
    fp->fgets(buf, 10000);
  } while(strncmp(buf, "***BEGIN PIPE OUTPUT:", 21));
  fscanf(fp->get_fp(), "%d\n", &nout);
  probs = new FLOAT[nout];
  lognout = (INT) ((FLOAT)logcode(1.0/(float)nout) * scale + 0.5);

  INT i;
  phone_threshold = new FLOAT[nout];
  if(floor_name == NULL) {
    for(i = 0; i < nout; i++)
      phone_threshold[i] = lna_floor;
  } else {
    pfp = new GENFILE(floor_name, "r");
    INT cnt = 0;
    for(i = 0; i < nout; i++){
      cnt += fscanf(pfp->get_fp(), "%f", &phone_threshold[i]);
      phone_threshold[i] = MIN(phone_threshold[i], lna_floor);
    }
    if(cnt != nout)
      panic("FILE_INTERFACE::initialize() failed to read %d phone floors from %s", 
	    nout, floor_name);
    delete pfp;
  }
    
  logpriors = new INT[nout];
  assert(nout == priors->size);
  for(i = 0; i < nout; i++)
    logpriors[i] = (INT) ((FLOAT)logcode(priors->get(i))*scale);

  if(Verbose) printf("Done (%d outputs)\n", nout);
}


VOID
BIN_FLOAT_INTERFACE::initialize(LIST<FLOAT> *priors, STRING floor_name, 
				STRING fname, FLOAT sc, 
				BOOL lf, INT nfl, INT ) {
  GENFILE *pfp;

  scale = sc;
  lower_floor = lf;
  n_for_lub_init = nfl;
  if(scale != 1.0){
    floor_prob = pow(prob_min, scale);
    floor_logcode = (INT) ((FLOAT)logcode_min * scale + 0.5);
  }  else {
    floor_prob = prob_min;
    floor_logcode = (INT)((FLOAT)logcode_min + 0.5);
  }
  if(Verbose) {
    printf("Initialising input probs (binary)...\n");
  }
  if(!strcmp(fname, "stdin") || !strcmp(fname, ""))
    fp = gstdin;
  else
    fp = new GENFILE(fname, "r");

  fp->fread((char*)&nout, sizeof(int), 1);
  probs = new FLOAT[nout];
  lognout = (INT) ((FLOAT)logcode(1.0/(float)nout) * scale + 0.5);

  INT i;
  phone_threshold = new FLOAT[nout];
  if(floor_name == NULL) {
    for(i = 0; i < nout; i++)
      phone_threshold[i] = lna_floor;
  } else {
    pfp = new GENFILE(floor_name, "r");
    INT cnt = 0;
    for(i = 0; i < nout; i++){
      cnt += fscanf(pfp->get_fp(), "%f", &phone_threshold[i]);
      phone_threshold[i] = MIN(phone_threshold[i], lna_floor);
    }
    if(cnt != nout)
      panic("FILE_INTERFACE::initialize() failed to read %d phone floors from %s", 
	    nout, floor_name);
    delete pfp;
  }
    
  logpriors = new INT[nout];
  assert(nout == priors->size);
  for(i = 0; i < nout; i++)
    logpriors[i] = (INT) ((FLOAT)logcode(priors->get(i))*scale);

  if(Verbose) printf("Done (%d outputs)\n", nout);
}

VOID
LNAFLOAT_INTERFACE::initialize(LIST<FLOAT> *priors, STRING floor_name, 
				STRING fname, FLOAT sc, 
				BOOL lf, INT nfl, INT no) {
  GENFILE *pfp;

  scale = sc;
  lower_floor = lf;
  n_for_lub_init = nfl;
  nout = no;
  if(scale != 1.0){
    floor_prob = pow(prob_min, scale);
    floor_logcode = (INT) ((FLOAT)logcode_min * scale + 0.5);
  }  else {
    floor_prob = prob_min;
    floor_logcode = (INT)((FLOAT)logcode_min + 0.5);
  }
  if(Verbose) {
    printf("Initialising input probs (binary)...\n");
  }

  if(!strcmp(fname, "stdin") || !strcmp(fname, ""))
    fp = gstdin;
  else
    fp = new GENFILE(fname, "r");

  //  fread((char*)&nout, sizeof(int), 1, fp);
  probs = new FLOAT[nout];
  lognout = (INT) ((FLOAT)logcode(1.0/(float)nout) * scale + 0.5);

  // lna floor and phone thresholds
  INT i;
  phone_threshold = new FLOAT[nout];
  if(floor_name == NULL) {
    for(i = 0; i < nout; i++)
      phone_threshold[i] = lna_floor;
  } else {
    pfp = new GENFILE(floor_name, "r");
    INT cnt = 0;
    for(i = 0; i < nout; i++){
      cnt += fscanf(pfp->get_fp(), "%f", &phone_threshold[i]);
      phone_threshold[i] = MIN(phone_threshold[i], lna_floor);
    }
    if(cnt != nout)
      panic("LNAFLOAT_INTERFACE::initialize() failed to read %d phone floors from %s", 
	    nout, floor_name);
    delete pfp;
  }
   
  // priors file 
  logpriors = new INT[nout];
  assert(nout == priors->size);
  for(i = 0; i < nout; i++)
    logpriors[i] = (INT) ((FLOAT)logcode(priors->get(i))*scale);

  if(Verbose) printf("Done (%d outputs)\n", nout);
}

VOID
LNA_INTERFACE::initialize(LIST<FLOAT> *priors, STRING floor_name, 
			  STRING fname, FLOAT sc, BOOL lf, 
			  INT nfl, INT no) {
  INT i;
  GENFILE *pfp;
  scale = sc;
  lower_floor = lf;
  n_for_lub_init = nfl;
  nout = no;
  if(scale != 1.0){
    floor_prob = pow(prob_min, scale);
    floor_logcode = (INT) ((FLOAT)logcode_min * scale + 0.5);
  }  else {
    floor_prob = prob_min;
    floor_logcode = (INT)((FLOAT)logcode_min + 0.5);
  }

  if(Verbose) {
    printf("Initialising input probs (lna)...\n");
  }

  if(!strcmp(fname, "stdin") || !strcmp(fname, "") || !strcmp(fname, "-"))
    fp = gstdin;
  else
    fp = new GENFILE(fname, "r");

  byte_buf = new UCHAR[nout+1];
  probs = new FLOAT[nout];
  prob_accum = new FLOAT[nout];
  data_priors = new FLOAT[nout];
  for(i = 0; i < nout; i++)
    prob_accum[i] = 0.0;
  lognout = (INT) ((FLOAT)logcode(1.0/(float)nout) * scale + 0.5);

  phone_threshold = new FLOAT[nout];
  if(floor_name == NULL) {
    for(i = 0; i < nout; i++)
      phone_threshold[i] = lna_floor;
  } else {
    pfp = new GENFILE(floor_name, "r");
    INT cnt = 0;
    for(i = 0; i < nout; i++){
      cnt += fscanf(pfp->get_fp(), "%f", &phone_threshold[i]);
      phone_threshold[i] = MIN(phone_threshold[i], lna_floor);
    }
    if(cnt != nout)
      panic("LNA_INTERFACE::initialize() failed to read %d phone floors from %s", 
	    nout, floor_name);
    delete pfp;
  }
    
  logpriors = new INT[nout];
  assert(nout == priors->size);
  for(i = 0; i < nout; i++) {
    data_priors[i] = priors->get(i);
    logpriors[i] = (INT) ((FLOAT)logcode(priors->get(i))*scale);
  }

  if(Verbose) printf("Done (%d outputs)\n", nout);
}

#ifdef SOCKETIO
VOID
SOCKET_INTERFACE::initialize(LIST<FLOAT> *priors, STRING floor_name, 
			     STRING ids, FLOAT sc, BOOL lf, 
			     INT nfl, INT ){
  INT n;
  INT portid = atoi(ids);
  GENFILE *pfp;

  scale = sc;
  lower_floor = lf;
  n_for_lub_init = nfl;
  if(scale != 1.0){
    floor_prob = pow(prob_min, scale);
    floor_logcode = (INT) ((FLOAT)logcode_min * scale + 0.5);
  }  else {
    floor_prob = prob_min;
    floor_logcode = (INT)((FLOAT)logcode_min + 0.5);
  }
  if(Verbose) 
    printf("Initialising probs_in socket()...\n");

  sock = new SOCKET();
  sock->open(portid);
  n = sock->read(&nout, sizeof(int));
  if(n == 0){
    panic("SOCKET_INTERFACE::initialize(%d): Premature end of probs in stream, reading nout\n", portid);
  } else if (n < 0){
    perror("SOCKET_INTERFACE::initialize()");
    panic("Read %d bytes for nout\n", -(n+1));
  } else {
    probs = new FLOAT[nout];
    lognout = (INT) ((FLOAT)logcode(1.0/(float)nout) * scale + 0.5);
    data_size = nout*sizeof(float);

    INT i;
    phone_threshold = new FLOAT[nout];
    if(floor_name == NULL) {
      for(i = 0; i < nout; i++)
	phone_threshold[i] = lna_floor;
    } else {
      pfp = new GENFILE(floor_name, "r");
      INT cnt = 0;
      for(i = 0; i < nout; i++){
	cnt += fscanf(pfp->get_fp(), "%f", &phone_threshold[i]);
	phone_threshold[i] = MIN(phone_threshold[i], lna_floor);
      }
      if(cnt != nout)
	panic("FILE_INTERFACE::initialize() failed to read %d phone floors from %s", 
	      nout, floor_name);
      delete pfp;
    }
    
    logpriors = new INT[nout];
    assert(nout == priors->size);
    for(i = 0; i < nout; i++)
      logpriors[i] = (INT) ((FLOAT)logcode(priors->get(i))*scale);

    if(Verbose) printf("Done (%d outputs)\n", nout);
  }
}
#endif

BOOL
ASCII_FLOAT_INTERFACE::read(INT *logprob, INT frame){
  INT id;
  BOOL skip = (logprob == NULL);

  INT nread = fscanf(fp->get_fp(), "%d\n", &id);
  if(id == END_OF_SENTENCE || fp->eof()) 
    return FALSE;
  if(nread == 0)
    panic("Unexpected end of input, frame = %d\n", frame);

  if(id != frame)
    panic("ASCII_FLOAT_INTERFACE::read(): frame mismatch, pipe %d, vs. dp %d\n",
	    id, frame);
  if(!skip)
    prob_sum = 0.0;
  for(INT i = 0; i < nout; i++){
    if(fscanf(fp->get_fp(), "%f", &probs[i]) != 1){
      perror("ASCII_FLOAT_INTERFACE::read(): Didn't read enough data.\n");
      panic("Read %d numbers, expecting %d\n", i, nout);
    }
    if(!skip) {
      prob_sum += probs[i];
      // convert from FLOAT probs to INT logprobs and subtract logprior
      if(lower_floor && probs[i] < phone_threshold[i]){
	probs[i] = floor_prob;
	logprob[i] = floor_logcode;
	nphones_pruned++;
      } else {
	logprob[i] = (INT) ((FLOAT)logcode(probs[i])*scale) - logpriors[i];
      }
    }
  }

  // set lub
  if(!skip) init_lub(frame);
  return(TRUE);
}

#ifdef ASCII_HEX
BOOL
ASCII_HEX_INTERFACE::read(INT *logprob, INT frame) {
  INT id;
  BOOL skip = (logprob == NULL);

  INT nread = fscanf(fp->get_fp(), "%d\n", &id);
  if(id == END_OF_SENTENCE || fp->eof()) 
    return FALSE;
  if(nread == 0)
    panic("Unexpected end of input, frame = %d\n", frame);

  if(id != frame)
    panic("ASCII_HEX_INTERFACE::read(): frame mismatch, pipe %d, vs. dp %d\n",
	    id, frame);

  if(!skip)
    prob_sum = 0.0;
  for(INT i = 0; i < nout; i++){
    if(fscanf(fp->get_fp(), "%x", &probs[i]) != 1){
      perror("ASCII_HEX_INTERFACE::read(): Didn't read enough data.\n");
      panic("Read %d numbers, expecting %d\n", i, nout);
    }
    if(!skip){
      prob_sum += probs[i];
      // convert from FLOAT probs to INT logprobs and subtract logprior
      if(lower_floor && probs[i] < phone_threshold[i]){
	probs[i] = floor_prob;
	logprob[i] = floor_logcode;
	nphones_pruned++;
      } else {
	logprob[i] = (INT) ((FLOAT)logcode(probs[i])*scale) - logpriors[i];
      }
    }
  }

  // set lub
  if(!skip) init_lub(frame);
  return(TRUE);
}
#endif

BOOL
BIN_FLOAT_INTERFACE::read(INT *logprob, INT frame) {
  INT id;
  BOOL skip = (logprob == NULL);

  INT nread = fp->fread((char*)&id, sizeof(int), 1);
  if(id == END_OF_SENTENCE || fp->eof()) 
    return FALSE;
  if(nread == 0)
    panic("Unexpected end of input, frame = %d\n", frame);

  if(id != frame)
    panic("BIN_FLOAT_INTERFACE::read(): frame mismatch, pipe %d, vs. dp %d\n",
	    id, frame);

  INT count;
  if((count = fp->fread((char*)probs, sizeof(FLOAT), nout)) != nout){
    perror("BIN_FLOAT_INTERFACE::read(): Didn't read enough data.\n");
    panic("Read %d numbers, expecting %d\n", count, nout);
  }

  if(!skip) {
    // convert from FLOAT probs to INT logprobs and subtract logprior
    prob_sum = 0.0;
    for(INT i = 0; i < nout; i++) {
      prob_sum += probs[i];
      if(lower_floor && probs[i] < phone_threshold[i]){
	probs[i] = floor_prob;
	logprob[i] = floor_logcode;
	nphones_pruned++;
      } else {
	logprob[i] = (INT) ((FLOAT)logcode(probs[i])*scale) - logpriors[i];
      }
    }
  
    // set lub
    init_lub(frame);
  }

  return(TRUE);
}

BOOL
LNAFLOAT_INTERFACE::read(INT *logprob, INT frame) {
  static INT id = 0;
  UCHAR label;
  BOOL skip = (logprob == NULL);

  if(id == END_OF_SENTENCE) {
    id = 0;
    return FALSE;
  }

  INT nread = fp->fread(&label, sizeof(UCHAR), 1);
  if(nread != 1) {
    if(fp->eof()) {
      id = END_OF_SENTENCE;
      return FALSE;
    } else {
      panic("LNAFLOAT_INTERFACE::read():Unexpected end of input, frame = %d\n", 
	    frame);
    }
  }
  if ((label & END_SENT_MASK) != 0)
    id = END_OF_SENTENCE;

  if((nread = fp->fread((char*)probs, sizeof(FLOAT), nout)) != nout){
    perror("LNAFLOAT_INTERFACE::read(): Didn't read enough data.\n");
    panic("Read %d numbers, expecting %d\n", nread, nout);
  }

  if(!skip) {
    // convert from FLOAT probs to INT logprobs and subtract logprior
    prob_sum = 0.0;
    for(INT i = 0; i < nout; i++){
      prob_sum += probs[i];
      if(lower_floor && probs[i] < phone_threshold[i]){
	probs[i] = floor_prob;
	logprob[i] = floor_logcode;
	nphones_pruned++;
      } else {
	logprob[i] = (INT) ((FLOAT)logcode(probs[i])*scale) - logpriors[i];
      }
    }

    // set lub
    if(frame >= 0)
      init_lub(frame);
  }

  return(TRUE);
}

BOOL
LNA_INTERFACE::read(INT *logprob, INT frame) {
  static INT id = 0;
  BOOL skip = (logprob == NULL);

  if(id == END_OF_SENTENCE) {
    id = 0;
    return FALSE;
  }

  if(lnaread(probs, byte_buf, nout, fp))
    id = END_OF_SENTENCE;

  if(fp->eof()) 
    return FALSE;

  if(!skip) {
    prob_sum = 0.0;
    for(INT i = 0; i < nout; i++){
      if(Debug_flag.phone)
	prob_accum[i] += probs[i];
      prob_sum += probs[i];
      // convert from FLOAT probs to INT logprobs and subtract logprior
      if(lower_floor && probs[i] < phone_threshold[i]){
	probs[i] = floor_prob;
	logprob[i] = floor_logcode;
	nphones_pruned++;
      } else {
	logprob[i] = (INT) ((FLOAT)logcode(probs[i])*scale) - logpriors[i];
      }
    }

    // set lub
    if(frame >= 0)
      init_lub(frame);
  }

  return(TRUE);
}

#ifdef SOCKETIO
BOOL
SOCKET_INTERFACE::read(INT *logprob, INT frame) {
  INT id = 0;
  UINT nread;
  BOOL skip = (logprob == NULL);

  nread = sock->read(&id, sizeof(int));

  if(nread == 0){
    id = END_OF_SENTENCE;
    sock->close();
  }

  /* return if end of sentence or end of data */
  if(id == END_OF_SENTENCE){
    return FALSE;
  }
  
  if(nread < sizeof(int)){
    sock->close();
    perror("SOCKET_INTERFACE::read()");
    panic("Unexpected end of input, frame = %d\n", frame);
  }
  
  if(id != frame){
    sock->close();
    panic("SOCKET_INTERFACE::read() mismatch frame %d (expecting %d)\n", id, frame);
  }

  /* read the data */
  nread = sock->read(probs, data_size);
  if(nread < data_size){
    sock->close();
    perror("Didn't read enough data from probs_in socket!");
    panic("Read %d bytes, expecting %d\n", nread, data_size);
  }

  if(!skip) {
    // convert from FLOAT probs to INT logprobs and subtract logprior
    prob_sum = 0.0;
    for(INT i = 0; i < nout; i++){
      prob_sum += probs[i];
      if(lower_floor && probs[i] < phone_threshold[i]){
	probs[i] = floor_prob;
	logprob[i] = floor_logcode;
	nphones_pruned++;
      } else {
	logprob[i] = (INT) ((FLOAT)logcode(probs[i])*scale) - logpriors[i];
      }
    }

    // set lub
    init_lub(frame);
  }

  return(TRUE);  
}
#endif

