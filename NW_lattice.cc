#include <ctype.h>
#include <math.h>
#include <string.h>

#include "NW_collection.h"
#include "NW_misc.h"
#include "NW_debug.h"
#include "NW_hypothesis.h"
#include "NW_lm.h"
#include "NW_lattice.h"

// DCA 06/MAR/96: order of operands changed: hypothesis counter added
LATTICE_LINK::LATTICE_LINK(INT t1, INT t2, INT w, INT a, INT v, INT r, INT l) {
  start_time = t1;
  stop_time = t2;
  start_node = stop_node = -1;
  wd = w;
  ac_prob = a;
  lm_prob = l;
  pron_ver = v; // DCA 06/MAR/96: variable name changed for transparency
  pron_prob = r;
  nlattice_elt++;
}

// DCA 06/MAR/96: changes required to print out both types of lattice
VOID 
LATTICE_LINK::print(GENFILE *fp, LIST<STRING> *vocab, INT index, INT nnodes){
  char* word;
  if(wd == Sentence_end_index)
    word = strdup(Silence_symbol);
  else
    word= strdup(vocab->get(wd & WORD_MASK));
    
  if(nnodes == -1){
    if(Use_pron_priors)
      fp->fprintf("J=%-5d  S=%-4d  E=%-4d W=%-19s v=%-2d a=%-12.6g r=%-12.6g\n", 
		index, start_node, stop_node, word, pron_ver, 
		logdecodelog(ac_prob), logdecodelog(pron_prob));
    else
      fp->fprintf("J=%-5d  S=%-4d  E=%-4d W=%-19s v=%-2d a=%-12.6g\n", 
		index, start_node, stop_node, word, pron_ver, 
		logdecodelog(ac_prob));
  } else {
    if(Use_pron_priors)
      fp->fprintf("J=%-5d S=%-4d E=%-4d W=%-19s v=%-2d a=%-12.6g  r=%-12.6g\n", 
		index, nnodes-start_node-1, nnodes-stop_node-1, word, pron_ver, 
		logdecodelog(ac_prob), logdecodelog(pron_prob)); 
    else
      fp->fprintf("J=%-5d S=%-4d E=%-4d W=%-19s v=%-2d a=%-12.6g\n", 
		index, nnodes-start_node-1, nnodes-stop_node-1, word, pron_ver, 
		logdecodelog(ac_prob)); 
  }
  delete word;
}

// DCA 06/MAR/96: changes required to print out both types of lattice
VOID 
LATTICE_LINK::print_with_frames(GENFILE *fp, LIST<STRING> *vocab, INT index){
  char* word;
  if(wd == Sentence_end_index)
    word = strdup(Silence_symbol);
  else
    word= strdup(vocab->get(wd & WORD_MASK));
    
    if(Use_pron_priors)
      fp->fprintf("J=%-5d  S=%-4d  E=%-4d W=%-19s v=%-2d a=%-12.6g r=%-12.6g\n", 
		  index, start_time, stop_time, word, pron_ver, 
		  logdecodelog(ac_prob), logdecodelog(pron_prob));
    else
      fp->fprintf("J=%-5d  S=%-4d  E=%-4d W=%-19s v=%-2d a=%-12.6g\n", 
		  index, start_time, stop_time, word, pron_ver, 
		  logdecodelog(ac_prob));

  delete word;

}

LATTICE_NODE::LATTICE_NODE(INT t, INT i) {
  time = t;
  if(i == -1) {
    id = nnodes++;
  } else {
    id = i;
    nnodes++;
  }
}

VOID 
LATTICE_NODE::print(GENFILE *fp, INT offset, FLOAT fshift, INT i) {
  if(i == -1) i = id;
  fp->fprintf("I=%d  t=%.3f\n", i, (FLOAT)(time+offset)*fshift);
}      
