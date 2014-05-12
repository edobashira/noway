// -*- Mode: C++;  -*-
// File: NW_linear_lexicon.h

class GRAMMAR;

#include "NW_acoustic_lexicon.h"

class PRONUNCIATION {
  QUEUE<SUBWORD*>* subwords;
public:
  PRONUNCIATION() {
    subwords = new QUEUE<SUBWORD*>();
  }

  PRONUNCIATION(QUEUE<SUBWORD*>* sq) { 
    subwords = sq->copy(); 
  }

  ~PRONUNCIATION() {
    delete subwords;
  }

  QUEUE<SUBWORD*>* get_subwords() { return subwords; }
};



class WORD {
  QUEUE<PRONUNCIATION*>* prons;
  INT word_id;
public:
  WORD(INT w){
    prons = new QUEUE<PRONUNCIATION*>();
    word_id = w;
  }

  ~WORD() {
    delete prons;
  }

  QUEUE<PRONUNCIATION*>* get_prons() { return prons; }
  INT get_word_id() { return word_id;}
  

  INT extend(INT t, HYP* h, QUEUE<STACK*> *stks);

};



class LINEAR_LEXICON : public ACOUSTIC_LEXICON {
  LIST<WORD*>* words;
  SUBWORD* interword_pause;
public:
  LINEAR_LEXICON() {
    words = new LIST<WORD*>();
  }

  ~LINEAR_LEXICON() {
    delete words;
  }

  LIST<WORD*>* get_words() { return words; }

  VOID build(PARAM_TABLE *parm, const INT nw, BOOL phi, GENFILE *lfp, 
    GENFILE *phfp, GENFILE *prfp = NULL);
  VOID score(HYP* h, GRAMMAR *lm, LIST<INT> *ps, LIST<INT> *ps2=NULL);
  INT extend(INT t, LIST<HYP*> *l, GRAMMAR *lm, QUEUE<STACK*> *stks);
  INT align(INT t, LIST<HYP*> *l, GRAMMAR *lm, QUEUE<STACK*> *stks);
};
