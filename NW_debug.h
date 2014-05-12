// -*- Mode: C++;  -*-
// File: NW_debug.h

// The global Debug_flag is used throughout to selectively
// print debug information.
class DEBUG_FLAG {
 public:
  DEBUG_FLAG() { clear(); }

  CHAR hyp;
  CHAR lub;
  CHAR bigram;
  CHAR trigram;
  CHAR cache;
  CHAR scheduler;
  CHAR text_decode;
  CHAR stack;
  CHAR extend;
  CHAR forward;
  CHAR prune;
  CHAR acoustic;
  CHAR phone;
  CHAR io;
  CHAR tree;
  CHAR tree_stats;
  CHAR decode_stats;
  CHAR how_are_we_doing;
  CHAR lm;
  CHAR phi;
  CHAR lattice;
  CHAR ndx;
  CHAR mixlm;

  VOID set_all() {
    hyp = lub = bigram = trigram = cache = prune = extend = scheduler = text_decode = stack = forward = acoustic = phone = io = tree = tree_stats = decode_stats = how_are_we_doing = lm =  phi = lattice = ndx = mixlm  = 1;
  }

  VOID clear() {
    hyp = lub = bigram = trigram = cache = extend = prune = scheduler = text_decode = stack = forward = acoustic = phone = io = tree = tree_stats = decode_stats = how_are_we_doing = lm = phi = lattice = ndx = mixlm = 0;
  }


#ifdef OLD_CPP
#define FLAG(NAME) if (!strcmp(flag,"NAME")) { NAME = 1; return; }
#define UNFLAG(NAME) if (!strcmp(flag,"NAME")) { NAME = 0; return; }
#else
#define FLAG(NAME) if (strcmp(flag, #NAME) == 0) { NAME = 1; return; }
#define UNFLAG(NAME) if (strcmp(flag, #NAME) == 0) { NAME = 0; return; }
#endif

  VOID set(CCHAR* flag) {
    FLAG(hyp)
    FLAG(extend)
    FLAG(forward)
    FLAG(prune)
    FLAG(bigram)
    FLAG(trigram)
    FLAG(cache)
    FLAG(scheduler)
    FLAG(text_decode)
    FLAG(stack)
    FLAG(phone)
    FLAG(acoustic)
    FLAG(lub)
    FLAG(io)
    FLAG(tree)
    FLAG(tree_stats)
    FLAG(decode_stats)
    FLAG(how_are_we_doing)
    FLAG(lm)
    FLAG(phi)
    FLAG(lattice)
    FLAG(ndx)
    FLAG(mixlm)
    printf("Unknown DEBUG_FLAG: %s\n", flag);
  }

  VOID unset(CCHAR* flag) {
    UNFLAG(hyp)
    UNFLAG(extend)
    UNFLAG(forward)
    UNFLAG(prune)
    UNFLAG(bigram)
    UNFLAG(trigram)
    UNFLAG(cache)
    UNFLAG(scheduler)
    UNFLAG(text_decode)
    UNFLAG(stack)
    UNFLAG(phone)
    UNFLAG(acoustic)
    UNFLAG(lub)
    UNFLAG(io)
    UNFLAG(tree)
    UNFLAG(tree_stats)
    UNFLAG(decode_stats)
    UNFLAG(how_are_we_doing)
    UNFLAG(lm)
    UNFLAG(phi)
    UNFLAG(lattice)
    UNFLAG(ndx)
    UNFLAG(mixlm)
    printf("Unknown DEBUG_FLAG: %s\n", flag);
  }

};

GLOBAL_VAR(DEBUG_FLAG, Debug_flag);
GLOBAL_VAR_INIT(GENFILE*, Decode_inf_file, NULL);
