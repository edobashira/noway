 // -*- Mode: C++;  -*-

#include "all.h"

// Phone model definition file format
// PHONE
// <nphones>
// <phone-id> <nstates> <phone-string>
// -1 -2 <distr-id> [nstates-2 times]
// <state> <ntrans-from-state> <to-state> <prob> ...
// ...
LIST<FLOAT> *
ACOUSTIC_LEXICON::make_phoneset(PARAM_TABLE *parm, BOOL phi, GENFILE *phfp, GENFILE *prfp) {
  CHAR buf[BUFSIZ];
  LIST<FLOAT> *priors = NULL;
  FLOAT dscale = parm->get_float("duration_scale", 1.0);
  FLOAT pdp = parm->get_float("phone_deletion_penalty", 1.0);

  if(phfp->fgets(buf, BUFSIZ) == NULL)
    panic("Failed to read from phone model file");
  if(phi == TRUE){
    if(strncmp(buf, "PHI", 3))
      panic("PHI file does not begin with PIF");
    priors = read_phi(parm, phfp, dscale, pdp);
  } else {
    if(strncmp(buf, "PHONE", 5))
      panic("Phone Models File does not begin with PHONE");
    read_phone_models(parm, phfp, dscale, pdp);
    priors = read_priors(prfp);
  }

  // build interword-pause (if not built)
  HMM* phone;
  if((phone = phoneset->get("interword-pause")) == NULL){
    FLOAT sprob, eprob;
    fprintf(stderr, "interword-pause not specified, automatically contructing...\n");
    if(silence_idx < 0) {
      if((phone = phoneset->get("sil")) != NULL){
	silence_idx = phone->distr;
      } else if ((phone = phoneset->get("h#")) != NULL){
	silence_idx = phone->distr;
      } else {
	silence_idx = 0;
	phone = phoneset->get(phonenames->get(0));
      }
    } else {
      phone = phoneset->get(phonenames->get(silence_idx));
    }
    fprintf(stderr, "Using phone %s for interword-pause\n", phone->name);
    sprob = phone->self_prob;
    eprob = phone->exit_prob;
    phone = new HMM(3, "interword-pause", silence_idx);
    phone->add_trans(0, 2, 1.0);
    phone->add_trans(phone->get_nstates()-1, phone->get_nstates()-1, sprob*dscale);
    if(parm->get_bool("pause_phone_deletion_penalty"))
      phone->add_trans(phone->get_nstates()-1, EXIT_STATE, eprob*dscale*pdp);
    else
      phone->add_trans(phone->get_nstates()-1, EXIT_STATE, eprob*dscale);
    phone->add_trans(0, 1, 1.0);

    // Store in hash table
    phonenames->push(strdup(phone->name));
    phoneset->insert(strdup(phone->name), phone);

    if(Debug_flag.phi)
      phone->print();
  } else {
    // phone = interword-pause
    if(phone->get_nstates() != 3)
      panic("interword-pause model should have ENTRY, EXIT and one non-null state!");
    if(!phone->check_trans(ENTRY_STATE, EXIT_STATE)){
      fprintf(stderr, "Adding (ENTRY_STATE, EXIT_STATE) transition to interword-pause model\n");
      phone->add_trans(ENTRY_STATE, EXIT_STATE, 1.0);
    }
  }

  return priors;
}

LIST<FLOAT> *
ACOUSTIC_LEXICON::read_phi(PARAM_TABLE *parm, GENFILE *phfp, 
			   FLOAT dscale, FLOAT pdp) {
  LIST<FLOAT> *priors = new LIST<FLOAT>;
  INT count= 0;
  INT nphones = 0;
  HMM* phone;
  CHAR buf[BUFSIZ];
  CHAR label[BUFSIZ];
  FLOAT pri, sprob, eprob;
  INT nstates;

  // PHI
  // # created on Mon Jul 15 06:08:42 1996 with:
  // # stateseq2phi -minFract 0.01 -ninp 15 -phntab /usr/local/abbot/etc/uklex/BEEPphone45.tab /home/spud5/ajr/PyneHillRobinson/PHR1.p16 tmp.phi
  // # number of phones:
  // 45
  // # label	prior		nstate	selfProb	exitProb
  // sil	0.469707	8	0.956890	0.043110
  // aa	        0.007041	4	0.868183	0.131817

  BOOL startflag = TRUE;
  while(phfp->fgets(buf, BUFSIZ) != NULL) {
    if(buf[0] == '#') {
      if(!strncmp(&buf[1], "SILENCE", 7))
	 sscanf(buf, "#SILENCE %d", &silence_idx);
      continue;
    }
    
    if(startflag) {
      startflag = FALSE;
      INT nph = sscanf(buf, "%d\n", &nphones);
      if(nph != 1) {
	fprintf(stderr, "Failed to read number of phones from PHI file");
      } else {
	continue;
      }
    }
    
    INT nr = sscanf(buf, "%s%f%d%f%f", label, &pri, &nstates, &sprob, &eprob);
    if(nr != 5)
      panic("failed to read PHI line:\n%s", buf);
    phone = new HMM(nstates+2, label, count);
    phone->phi_init(sprob, eprob, dscale, pdp);
    if(Debug_flag.phi)
      phone->print();

    // Store in hash table
    phonenames->push(strdup(label));
    phoneset->insert(strdup(label), phone);

    priors->push(pri);

    count++;
  }
  if(nphones == 0)
    nphones = count;
  else
    assert(nphones == count);

  return priors;
}

VOID
HMM::phi_init(FLOAT sprob, FLOAT eprob, FLOAT dur_scale, FLOAT pdp){
  self_prob = sprob;
  exit_prob = eprob;

  // initial transition
  add_trans(ENTRY_STATE, FIRST_STATE, 1.0);

  // non-final real states
  for(int i = FIRST_STATE; i < nstates-1; i++){
    if(Forward_process) {
      // Forward prob model has 0.5/0.5 transitions
      add_trans(i, i+1, 0.5*dur_scale);
      add_trans(i, i, 0.5*dur_scale);
    } else {
      // Viterbi model has 1.0 probs until final state
      add_trans(i, i+1, 1.0*dur_scale);
    }
  }

  // final state selfloop
  add_trans(nstates-1, nstates-1, self_prob*dur_scale);

  // exit transition
  add_trans(nstates-1, EXIT_STATE, exit_prob*dur_scale*pdp);
}

LIST<FLOAT> *
ACOUSTIC_LEXICON::read_priors(GENFILE *prfp) {
  LIST<FLOAT> *priors = new LIST<FLOAT>;
  FLOAT p;
  while(fscanf(prfp->get_fp(), "%f", &p) == 1)
      priors->push(p);
  if(!prfp->eof())
    panic("Premature end of priors file after %d priors", priors->size);
  return priors;
}
  
VOID
ACOUSTIC_LEXICON::read_phone_models(PARAM_TABLE *parm, GENFILE *phfp, 
				    FLOAT dscale, FLOAT pdp) {
  UINT nphones;

  INT nread = fscanf(phfp->get_fp(), "%d", &nphones);
  if(nread != 1)
    panic("Failed to read in num. of phones from Phone Model File");
  
  //  Read in the phone model defs, build the model objects and
  //   store in a hash table
  INT state, state2, nstates, id, j, k, ntrans, dis;
  FLOAT prob;
  CHAR name[BUFSIZ];
  HMM* phone;
  for(UINT i = 0; i < nphones; i++){
    // Read in phone name and number of states
    nread = fscanf(phfp->get_fp(), "%d%d%s", &id, &nstates, name);
    if(nread != 3)
      panic("Error reading phone information (%s (%d), %d states)\n",
	    name, id, nstates);
    nread = fscanf(phfp->get_fp(), "%d%d", &state, &state2);
    if (nread != 2)
      panic("Error in reading %s phone state IDs", name);
    if(state != -1 || state2 != -2)
      panic("Expecting -1(%d) and -2(%d) for null states in %s", 
	    state, state2, name);
    for(j = 0; j < nstates-2; j++) {
      nread = fscanf(phfp->get_fp(), "%d", &dis);
      if(nread != 1)
	panic("Unexpected end of %s phone id input", name);
    }

    // Create new phone object
    phone = new HMM(nstates, name, dis);

    // Write in transitions
    for(j = 0; j < nstates; j++){
      nread = fscanf(phfp->get_fp(), "%d", &state);
      if(nread != 1)
	panic("Phone %s:  Failed reading transitions (state %d)\n", name, j);
      if(state != j)
	panic("Phone %s:  Error in reading transitions (expecting s %d, got %d)\n",
	      name, j, state);
      nread = fscanf(phfp->get_fp(), "%d", &ntrans);
      if(nread != 1)
	panic("Phone %s: Failed reading ntrans from state %d\n", name, state);
      for(k = 0; k < ntrans; k++){
	nread = fscanf(phfp->get_fp(), "%d%f", &state2, &prob);
	if(nread != 2)
	  panic("Phone %s: Failed reading trans from state %d", name, state);
	// phone deletion penalty for transitions to end state
	if(state2 == EXIT_STATE && state != ENTRY_STATE) {
	  // only want phone deletion penalty on transitions from 
	  // non-null states to EXIT_STATE
	  if(parm->get_bool("pause_phone_deletion_penalty") ||
	     strcmp(name, "interword-pause"))
	    // only want phone deletion penalty in interword-pause
	    // if -pause_phone_deletion_penalty specified
	    prob = prob * pdp;
	}
	// Duration scaling if required
	phone->add_trans(state, state2, prob*dscale);
      }
      //      if(phone->trans_prob == NULL || phone->ntrans == 0)
      if(phone->trans == NULL || phone->ntrans == 0)
	panic("No transitions given for phone %s\n", name);
      if(phone->get_trans()[0].from != ENTRY_STATE)
	panic("No initial transition for phone %s\n", name);
    }

    // Store in hash table
    phonenames->push(strdup(name));
    phoneset->insert(strdup(name), phone);

  }



  nread = fscanf(phfp->get_fp(), "%d", &id);
  if(nread != EOF)
    panic("Failed to reach end of file after %d phones", nphones);
  assert(phoneset->size() == nphones);
}



VOID
ACOUSTIC_LEXICON::set_acoustic_input(PARAM_TABLE *p, LIST<FLOAT> *priors){

#ifdef ASCII_HEX
#ifdef SOCKETIO
  assert(p->get_bool("lna") || p->get_bool("lnafloat") || p->get_bool("ascii") ||
	 p->get_bool("bin") || p->get_bool("hex") || p->get_bool("socket"));
#else
  assert(p->get_bool("lna") || p->get_bool("lnafloat") || p->get_bool("ascii") ||
	 p->get_bool("bin") || p->get_bool("hex"));
#endif
#else
#ifdef SOCKETIO
  assert(p->get_bool("lna") || p->get_bool("lnafloat") || p->get_bool("ascii") ||
	 p->get_bool("bin") || p->get_bool("socket"));
#else
  assert(p->get_bool("lna") || p->get_bool("lnafloat") || p->get_bool("ascii") ||
	 p->get_bool("bin"));
#endif
#endif
  STRING pfloorf = NULL;
  if(p->get_bool("phone_floor"))
    pfloorf = strdup(p->get_string("phone_floor", "phone.floors"));
  FLOAT as_default = 1.0;
  if(p->get_bool("demo")) {
    as_default = 0.3;
    if(Verbose && ! p->get_bool("acoustic_scale"))
      printf("-demo option setting -acoustic_scale 0.3\n");
  }
  FLOAT scale = p->get_float("acoustic_scale", as_default);
  INT nfl = p->get_int("lub_n", 5);

  pause_model->prob_min = 1.0e-10;
  BOOL lf = p->get_bool("prob_min") || p->get_bool("phone_floor") || p->get_bool("demo");
  if(lf){
    if(p->get_bool("prob_min") || p->get_bool("demo")) {
      pause_model->lna_floor = p->get_float("prob_min", 0.00025, 0);
      pause_model->prob_min = p->get_float("prob_min", 1.0e-10, 1);
      if(Verbose && p->get_bool("demo"))
	printf("-demo option setting -prob_min 0.00025\n");
    } else {
      pause_model->lna_floor = 1.0;
    }
  } else {
    pause_model->lna_floor = 0.0;
  }
  pause_model->logcode_min = logcode(pause_model->prob_min);

  //  STRING priorf = strdup(p->get_string("priors", "priors.WSJ"));
  INT nou;
  if(p->get_bool("lna")) {
    if((nou = p->get_int("lna", priors->size, 0)) != priors->size)
      panic("Specified -lna %d expecting -lna %d\n", nou, priors->size);
    pause_model->lna_in(priors, pfloorf, p->get_string("lna", "stdin", 1), 
		 scale, lf, nfl, nou);
  } else if (p->get_bool("lnafloat")) {
    if((nou = p->get_int("lnafloat", priors->size, 0)) != priors->size)
      panic("Specified -lnafloat %d expecting -lnafloat %d\n", nou, priors->size);
    pause_model->lnafloat_in(priors, pfloorf, p->get_string("lnafloat", "stdin", 1), 
		      scale, lf, nfl, nou);
  } else if (p->get_bool("ascii"))
    pause_model->ascii_float_in(priors, pfloorf, p->get_string("ascii", "stdin"), scale, lf, nfl);
  else if (p->get_bool("bin"))
    pause_model->bin_float_in(priors, pfloorf, p->get_string("bin", "stdin"), scale, lf, nfl);
#ifdef ASCII_HEX
  else if (p->get_bool("hex"))
    pause_model->ascii_hex_in(priors, pfloorf, p->get_string("hex", "stdin"), scale, lf, nfl);
#endif
#ifdef SOCKETIO
  else if (p->get_bool("socket"))
    pause_model->socket_in(priors, pfloorf, p->get_int("socket", 5765), scale, lf, nfl);
#endif
  else
    panic("ACOUSTIC_LEXICON::set_acoustic_input():failed to specify input stream");

  if(p->get_bool("no_lna_check"))
    pause_model->unset_sum_check();
  else
    pause_model->set_sum_check();
}


