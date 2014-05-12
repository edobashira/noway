// -*- Mode: C++;  -*-

class LATTICE_LINK {
private:
//  INT link_id; // DCA 07/MAR/96: variable name changed for transparency
  INT start_time;
  INT stop_time;
  INT start_node;
  INT stop_node;
  INT wd;
  INT pron_ver; // DCA 06/MAR/96: variable name changed for transparency
  INT ac_prob;
  INT lm_prob;
  INT pron_prob;
public:
  LATTICE_LINK(INT t1, INT t2, INT w, INT a, INT v = 0, INT r = 0, INT l = 0);
  ~LATTICE_LINK() { if(!pipe) nlattice_elt--; }

  VOID reset() { nlattice_elt = 0; }
  VOID inc() {nlattice_elt++;}

  INT get_start_time() { return(start_time); }
  INT get_end_time() { return(stop_time); }
  INT get_start_node() { return(start_node); }
  INT get_end_node() { return(stop_node); }
  INT get_word() { return wd; }
  INT get_version() { return pron_ver; }
  INT get_aprob() { return ac_prob; }
  INT get_lprob() { return lm_prob; }
  INT get_pron_prob() { return pron_prob; }

  VOID put_start_node(INT nid) { start_node = nid; }
  VOID put_end_node(INT nid) { stop_node = nid; }

  static INT nlattice_elt;
  static BOOL pipe;
  INT num_links() { return(nlattice_elt); }


  VOID print(GENFILE *fp, LIST<STRING> *vocab, INT index, INT nnodes = -1);
  VOID print_with_frames(GENFILE *fp, LIST<STRING> *vocab, INT index);
};  


class LATTICE_NODE {
private:
  INT id;
  INT time;
  static INT nnodes;
public:

  LATTICE_NODE(INT t, INT i = -1);
  ~LATTICE_NODE() { nnodes--; }

  VOID print(GENFILE *fp, INT offset = 0, FLOAT fshift=0.016, INT i = -1);
  
};

