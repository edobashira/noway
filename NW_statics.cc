// -*- Mode: C++;  -*-
// File: statics.cc
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
//* Created: Wed Apr  6 14:47:36 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "all.h"
#include "NW_scheduler.h"

INT          HMM::ndistr = 0;
INT          HMM::nmodels = 0;
INT          HMM::total_states = 0;
QUEUE_OFFSET<INT*>*  ACOUSTIC::outputs = new QUEUE_OFFSET<INT*>(2048);
QUEUE_OFFSET<INT>*   ACOUSTIC::lub     = new QUEUE_OFFSET<INT>(2048);
QUEUE_OFFSET<INT>*   ACOUSTIC::def_lub_inc = new QUEUE_OFFSET<INT>(2048);
INT          ACOUSTIC::sentence_length = INT_MAX;
INT          ACOUSTIC::sentence_count = -1;
INT          ACOUSTIC::frame_count = 0;
INT          ACOUSTIC::frames_this_sentence = 0;
INT          ACOUSTIC::proc_frame_count = 0;
INT          ACOUSTIC::proc_frames_this_sentence = 0;
INT          ACOUSTIC::nout = 0;
INT          ACOUSTIC::logtrans = logcode(0.5);
INT          ACOUSTIC::beam = 0;
INT          ACOUSTIC::state_beam = 0;
INT          ACOUSTIC::pathword_threshold = 0;
INT          ACOUSTIC::start_time = 0;
INT          ACOUSTIC::nnodes_used = 0;
INT          ACOUSTIC::nphones_pruned = 0;
//INT          ACOUSTIC::nreactivated = 0;
FLOAT        ACOUSTIC::lna_floor = 0.000025;
FLOAT        ACOUSTIC::prob_min = 1.0e-10;
INT          ACOUSTIC::logcode_min = -655360;
BOOL         ACOUSTIC::strict_lub = FALSE;
INT          ACOUSTIC::startprob = 0;
INT          ACOUSTIC::garbage_model = -1;
INT          ACOUSTIC::average_lm_contexts = FALSE;
INT          ACOUSTIC::lmscale = 1;
INT          ACOUSTIC::lmoffset= 0;
FLOAT        ACOUSTIC::amscale = 1.0;
FLOAT        ACOUSTIC::logbase = 10.0;
FLOAT*       ACOUSTIC::prob_accum = NULL; 
FLOAT*       ACOUSTIC::data_priors = NULL; 
BOOL         ACOUSTIC::ndx_output = FALSE;
INT          LATTICE_LINK::nlattice_elt = 0;
INT          LATTICE_LINK::pipe = FALSE;
INT          LATTICE_NODE::nnodes = 0;
ACOUSTIC_INTERFACE* HMM::input = NULL;
INT          NODE::nnodes = 0;
BOOL         NODE::smear_unigram = FALSE;
BOOL         NODE::smear_context = TRUE;
QUEUE<INT>*  NODE::offp_queue = new QUEUE<INT>(128);
LIST<SUBWORD*>* SUBWORD::active_sw_list = new LIST<SUBWORD*>(1024);
INT          SUBWORD::nsubwords = 0;
BACKTRACE_LIST* BACKTRACE_LIST::free_list = NULL;
INT          BACKTRACE_LIST::block_size = 1024;
INT          BACKTRACE_LIST::total_backtrace_lists = 0;
  
INT          STACK::max_hyps = 0;
BOOL         STACK::merge_mode = FALSE;
//INT          TRIGRAM_OBJ::ncreated = 0;
//INT          BIGRAM_OBJ::ncreated = 0;
//INT          UNIGRAM_OBJ::ncreated = 0;
HYP*         HYP::free_list = NULL;
INT          HYP::block_size = 2400;
INT          HYP::total_hyps = 0;
INT          HYP::lm_context_length = 3;
INT          HYP::mixture_lm_size = 0;
BOOL         HYP::output_pron_version = FALSE;
BOOL         HYP::output_silence = FALSE;
INT          HYP::lm_prob = 0;
LIST<INT>*   HYP::gamma = NULL;
LIST<INT>*   HYP::ng_probs = NULL;

// parameter options
STRING       SCHEDULER::decoder_param[] = {
  "dummy command name",  // a hack to allow this to be read as a argc/argv style command line
  "-parameter_files", "list of file names to read as parameters",
// debug and vebose options
  "-version", "print version information and quit",
  "-verbose", "print occasional diagnostics",
  "-help", "print out list of all paramerers",
  "-d", "debug options",
  "-n_hyps_db", "number of hypotheses to display each time with -d stack (5)",
  "-max_sents", "only decode this many sentences",
  "-skip_sents", "skip the first <n> sentences",
  "-skip_frames", "skip the first <s> frames in the first sentence to be decoded",
  "-online", "online processing",
  "-reset_time_each_utt", "timings for each utt start from 0",
// Grammar stuff
  "-blind_mixture_ngram", "mixture of ngrams file",
  "-lsa_mixture_ngram", "lsa mixture LM file",
  "-ngram", "n-gram grammar file",
  "-trigram", "trigram grammar file",
  "-bigram", "bigram grammar file",
  "-tagprobs", "unigram tag probabilities", 
  //  "-wordpair", "wordpair grammar file",
  "-nogram", "don't use a grammar",
  "-fsn", "finite state network grammar",
  "-lm_scale", "language model scale value (multiplies LM logprobs)",
  "-lm_offset", "language model offset value (added to LM logprobs)",
  "-use_unk", "use <UNK> in LM",
  "-unk_factor", "factor to scale P(<UNK> | *)",
  //  "-large_vocab", ">64K word dictionary", // not yet done
  "-write_lm", "write out a binary LM file",
  //  "-byte_swap_write_lm", "byte swap when writing binary LM file",
  "-print_lm", "print out LM file",
  "-print_unigrams", "print out LM unigrams",
  "-cross_sentence", "cross sentence decoding (i.e. ... foo </s> <s> bar ...)",
  //  "-lm_merge_param", "smoothing parameter for different LM estimate",
  //  "-lm_merge_log_domain", "merge in log LM probability domain",
  //  "-unigram_merge", "reestimated unigram probs",
  "-average_lm_contexts", "Average LM probs for all hyps",
// Phone model stuff
  "-phone_models", "phone models file",
  "-priors", "phone prior probability file",
  "-phi", "phone information file (this may be unreliable; see the man page)",
  "-duration_scale", "duration scale factor",
  "-phone_deletion_penalty", "multiplicative scaling of phone exit probs",
  "-pause_phone_deletion_penalty", "use phone deletion penalty for interword pause",
  "-garbage_model", "Garbage model 'phone'",
// Dictionary stuff
  "-dictionary", "pronunciation dictionary file",
  "-no_pron_prior", "ignore priors on multiple pronunciations", 
  "-missing_pron_file", "missing pronunciations go here", 
  "-sentence_start_symbol", "symbol for beginning of sentence (<s>)",
  "-sentence_end_symbol", "symbol for end of sentence (</s>)",
  "-paragraph_symbol", "symbol for paragraph start (<p>)",
  "-article_symbol", "symbol for article start (<art>)",
  "-silence_symbol", "symbol for article start (<art>)",
  "-unknown_word_symbol", "symbol for words unkown to LM (<UNK>)",
// Acoustic probabilities stuff
  "-lna", "lna-coded phone probability file(s)",
  "-lnafloat", "32 bit float phone probability file",
  "-ascii", "ascii float phone probability file",
  "-bin", "binary float phone probability file",
#ifdef ASCII_HEX
  "-hex", "ascii hex phone probability file",
#endif
#ifdef SOCKETIO
  "-socket", "port id of socket communicating phone probabilites",
#endif
  "-acoustic_scale", "language model acoustic match factor",
  "-no_lna_check", "don't check that lna files sum to 1",
  "-frame_shift", "frame shift of acoustic analysis in ms (16)",
// Pruning stuff
  "-demo", "Sets beam/state_beam/prob_min/n_hyps/smear_unigram/new_lub",
  "-prob_min", "floor value for acoustic posteriors (<floor> <replacement>)",
  "-phone_floor", "phone dependent floor values",
  "-n_hyps", "number of active hypotheses to consider starting at each time",
  "-pathword_threshold", "max num of pathwords in node to consider LM probs",
  "-beam", "beam around lub, outside of which prune (word ends)",
  "-state_beam", "beam around lub, outside of which prune (states)",
  "-smear_unigram", "smear unigram LM probs in tree",
  "-smear_context", "do not use max_w p(.|.,w) probabilities as default LM",
  "-hyp_prune", "prune hypotheses in the tree", 
  "-lub_n", "sum over probabilities 1...n to initialise lub(t)",
  "-new_lub", "only allow backtraces from word ends to add to lub",
// Decoding output
  "-output_file", "file to output decoding in [stdout]",
  "-inc_output", "write incremental output to stdout",
  "-utterance_score", "print out logprob following decoding",
  "-lattice", "write out a word lattice (HTK format)",
  //  "-detailed_lattice", "write out a more detailed lattice (HTK format)",
  "-pipe_lattice", "output lattice elements as they are created", 
  "-ctm", "output in ctm format", 
  "-srt", "output in NIST / TREC Speech Retrieval format",
  "-output_pron_version", "output multiple pronunciation information in CTM/SRT/log files",
  "-output_silence", "output silence segments in CTM/SRT files",
  "-ndx", "Index file for Broadcast News decoding (should correspond to lna)",
  "-utt_desc_file", "utterance description file (1 line per utterance)",
  "-phone_decode", "state path decoding",
  "-dump_lub", "write out a file for plotting of lub and best hypothesis",
  //  "-forced_alignment", "perform forced alignments",  
// Decoding algorithm
  "-align", "perform forced alignment",
  "-merge_hyps", "merge hypotheses at word ends",
  "-forward_process", "do forward processing (not Viterbi) within a word",
  "-no_decode", "no decoding performed (used when reading/writing LM only)",
  "-text_decode", "LM score some text", 
  "-nl_is_end_sent", "Treat newline as </s><s> when text-decoding",
// Optimizations
  "-no_lm_cache", "don't cache language model",
  //  "-full_lm_cache", "full cache of language model",
  //  "-old_caching", "Use old incremental caching scheme",
// End of param table
  NULL, NULL,
};


