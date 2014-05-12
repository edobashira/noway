// -*- Mode: C++;  -*-
// File: ngram.cc
// Author: Steve Renals (sjr@dcs.shef.ac.uk)
// Copyright (C) Department of Computer Science, Sheffield University, 1997
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* REQUIRES:
//*
//* HISTORY:
//* Created: 
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_debug.h"
#include "NW_param.h"

#include "NW_hypothesis.h"
#include "NW_fsn.h"
#include "NW_stack.h"
#include "NW_lexicon.h"

INT 
FSN::build(GENFILE *, LEXICON *, INT , INT , BOOL , INT , BOOL){
  panic("FSN::build() not yet implemented!");
  return(0);
}

// Builds from a linear ref
INT 
FSN::build_from_ref(QUEUE<INT>* wq){
  FSN_NODE* start = root;
  FSN_NODE* end;
  INT n = 0;
  while(wq->size > 0){
    n++;
    end = new FSN_NODE();
    start->get_links()->push(new FSN_LINK(wq->dequeue(), start, end));
    end->inc_nin();
    start = end;
  }
  return n;    
}


LIST<INT>*
FSN::get_successors(HYP* h){
  LIST<INT>* res = new LIST<INT>();
  FSN_NODE* nd = current_node(h);
  if(nd == NULL)
    return res;
  LIST<FSN_LINK*>* l = nd->get_links();
  l->init_cursor();
  while(l->next()) 
    res->push(l->read()->get_word());
  return res;
}


INT
ALIGN_LM::build(GENFILE *fp, LEXICON *lex, INT, INT, BOOL, INT, BOOL) {
  QUEUE<INT>* wq = new QUEUE<INT>();
  INT sent = 0;
  INT w;

  while(fp->fgets(albuffer, ALIGN_BUFFER_SIZE) != NULL) {
    sent++;
    if(strlen(albuffer) == (UINT)ALIGN_BUFFER_SIZE-1)
      panic("Buffer size exceeded! Max %d characters in an utterance ref\n",
	    ALIGN_BUFFER_SIZE-1);

    wq->fast_clear();
    CHAR* refptr = strtok(albuffer, " \t");
    do {
      if(refptr[0] == '(' && refptr[strlen(refptr)-1] == ')')
	break;  // utterance id
      w = lex->vocab_hash->get(refptr);
      if(w == -1){
	w = lex->vocabulary->size;
	lex->vocab_hash->insert(strdup(refptr), w);
	lex->vocabulary->push(strdup(refptr));
      }
      if(!strcmp(refptr, Silence_symbol))
	Silence_index = w;
      wq->enqueue(w);
    } while((refptr = strtok(NULL, " \t")) != NULL);

    if(refptr != NULL){
      refptr = strtok(NULL, " \t");
      if(refptr != NULL)
	panic("More text for ref %d after sentence id!\n", sent);
    }
    
    FSN* fsn = new FSN();
    fsn->build_from_ref(wq);
    utt_models->enqueue(fsn);
  }
  return lex->vocabulary->size;
}
