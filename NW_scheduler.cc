// -*- Mode: C++;  -*-
// File: scheduler.cc
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
// Copyright (C) Sheffield University, 1995
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* REQUIRES:
//*
// -*- Mode: C++;  -*-
// File: scheduler.cc
// Author: Steve Renals (sjr@dcs.shef.ac.uk)
// Copyright (C) Department of Computer Science, Sheffield University, 1995-98
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* REQUIRES:
//*
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* HISTORY:
//*  Apr 20 10:38 1995 (sjr): Empty stack printout condition caught  
//*                            (caused seg fault with narrow beam)
//*  Jan 27 13:00 1995 (sjr): Supports state path decoding
//* Created: Fri Apr 22 13:30:39 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#include "all.h"
#include "NW_scheduler.h"
#ifdef TIMING
#include <sys/utsname.h>
#endif

VOID 
SCHEDULER::main_prog(INT argc, STRING argv[]){
  STACK *transcriptions = new STACK();

  welcome_mesg();

  param = new PARAM_TABLE(NULL, decoder_param);
  parse(argc, argv);
  if(Verbose){
#ifdef TIMING
    time(&tt);
    CHAR hname[BUFSIZ];
    struct utsname unamestruct;
    uname(&unamestruct);
    strcpy(hname, unamestruct.nodename);
    printf("%s: %s\n", hname, ctime(&tt));
#endif
    param->print();
    fflush(stdout);
  }
  sentence_count = 1;
  maxhyps = param->get_int("n_hyps", 7);
  if(Verbose && param->get_bool("demo"))
    printf("-demo option setting -n_hyps %d\n", maxhyps);
  build();
  transcriptions->set_max_hyps(maxhyps);
  if(param->get_bool("merge_hyps"))
    transcriptions->set_merge_mode();
  n_hyps_db = param->get_int("n_hyps_db", 5);
  output_logprob = param->get_bool("utterance_score");
#ifdef TIMING
  endt = times(&end_time_buffer);
  if(Verbose){
    printf("Setup Timing Information:\n");
    printf("  Elapsed real time = %lds\n", (endt-startt)/tick);
    printf("  User time = %lds\n", 
	   (end_time_buffer.tms_utime-start_time_buffer.tms_utime)/tick);
    printf("  System time = %lds\n", 
	   (end_time_buffer.tms_stime-start_time_buffer.tms_stime)/tick);
    fflush(stdout);
  }
#endif
  if(Debug_flag.decode_stats) {
    Decode_inf_file = new GENFILE("DECODE.STATS", "w");
  }
  
#ifdef TIMING
  startt = times(&start_time_buffer);
#endif
  delete transcriptions;

  INT max_sent = INT_MAX;
  if(param->get_bool("max_sents"))
    max_sent = param->get_int("max_sents") + sentence_count - 1;
    
  while((transcriptions = decode()) != NULL){
    output_hyps(transcriptions);
    delete transcriptions;
    if(++sentence_count > max_sent) break;
  }

  if(ndx_output){
    if(srt_output){
      srt_fp->fprintf("%s\n", ndx_ftr);
      delete srt_fp;
    }
    if(ctm_output){
      ctm_fp->fprintf("; %s\n", ndx_ftr);
      delete ctm_fp;
    }
  }
    
#ifdef TIMING
  sentence_tick = (sentence_count-1)*tick;
  endt = times(&end_time_buffer);
#endif
  if(Verbose){
    printf("\nDecoding Statistics:\n");
    if(Debug_flag.phone && prob_accum != NULL) {
      printf("Average posteriors vs. priors:\n");
      for(INT i = 0; i < nout; i++)
	printf(" %s: %f   %f   (%f)\n", dictionary->phonenames->get(i), 
	       prob_accum[i]/(FLOAT)total_proc_frames_seen(), data_priors[i],
	       prob_accum[i]/((FLOAT)total_proc_frames_seen() * data_priors[i]));
    }
#ifdef TIMING
    printf("  Elapsed real time = %lds (%lds/sentence)\n", 
	   (endt-startt)/tick, (endt-startt)/sentence_tick);
    printf("  User time = %lds (%lds/sentence)\n", 
	   (end_time_buffer.tms_utime-start_time_buffer.tms_utime)/tick,
	   (end_time_buffer.tms_utime-start_time_buffer.tms_utime)/sentence_tick);
    printf("  System time = %lds (%lds/sentence)\n", 
	   (end_time_buffer.tms_stime-start_time_buffer.tms_stime)/tick,
	   (end_time_buffer.tms_stime-start_time_buffer.tms_stime)/sentence_tick);
    printf("  User Time Speed = %.2f xrealtime\n", 
	   (FLOAT) ((end_time_buffer.tms_utime-start_time_buffer.tms_utime)/tick)
	   / (frame_shift*total_proc_frames_seen()));
    printf("  Wallclock Speed = %.2f xrealtime\n\n", 
	   (FLOAT) ((endt-startt)/tick) / (frame_shift*total_proc_frames_seen()));
#endif
    if(Debug_flag.tree_stats) {
      printf("  Histogram of node activations by tree depth:\n");
      INT depth = 0;
      while(depth_count->size >0){
	printf("    Level %d: %d nodes activated\n", depth, depth_count->dequeue());
	depth++;
      }
    }
    printf("  Average active models/frame: %ld\n", total_nodes_used/total_proc_frames_seen());
    printf("  Maximum active models/start: %d\n", max_nnodes_used);
    if(Debug_flag.tree_stats)
      printf("    (Average %ld above pw threshold)\n", total_nodes_pw/total_proc_frames_seen());
    printf("  Average unique models/start: %ld\n", total_nodes_used2/total_proc_frames_seen());
    printf("  Maximum unique models/start: %d\n", max_nnodes_used2);
    printf("  Average number of hypotheses created/start: %ld\n",
	   total_hyps_created/total_proc_frames_seen());
    printf("  Maximum number of hypotheses created/start: %d\n",
	   max_nhyps_created);
    printf("  Average number of hypotheses used/start: %ld\n",
	   total_hyps_used/total_proc_frames_seen());
    printf("  Average activation duration: %.1f\n", 
	   (FLOAT)total_nodes_used/(FLOAT)total_nodes_used2);
    printf("  Average phones pruned by PDP: %.1f%%\n", 
	   (FLOAT)(100.0 * total_phones_pruned) / (FLOAT)(total_proc_frames_seen() * nout));
    fflush(stdout);
  }

#if 0
#if defined(Linux) || defined(Solaris2)
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    printf("Memory Usage:\n");
    printf("  Max RSS = %ld, I RSS = %ld,  Data = %ld, Stack = %ld\n",
	   ru.ru_maxrss, ru.ru_ixrss, ru.ru_idrss, ru.ru_isrss);
    printf("  Swaps = %ld,  Page reclaims = %ld,  Page faults = %ld\n",
	   ru.ru_nswap, ru.ru_minflt, ru.ru_majflt);
#endif /*Solaris2 Linux*/
#endif /* 0 */

  if(Debug_flag.decode_stats) {
    delete Decode_inf_file;
  }
}

// Returns n-best hypotheses in a priority queue;  returns NULL if end
//  of data
STACK*
SCHEDULER::decode() {
  LIST<HYP*>* hlist = new LIST<HYP*>();
  LIST<HYP*>* delhlist = new LIST<HYP*>();
  STACK *pq, *res;
  HYP *h = NULL;
  INT nh;
  INT top_prob;
  INT nhyps_created, nhyps_used;
  INT nnodes_used_sum, nnodes_used_sum2, nnodes_pw;
  INT dl;

  INT save_beam = beam;
  INT save_sbeam = state_beam;
//  BOOL reset_outputs;
//  BOOL restart, restarted;
  UINT smask = 0x10000;
  UINT sunmask = 0xffff;
  UINT maxmask = 0x80000000;

  if(ndx_output) {
    if(ndxq->size == 0)
      return (STACK*)NULL;
    ndx_sec = ndxq->dequeue();
    INT nskip = ndx_sec->get_start_time() - total_frames_seen();
    if(nskip < 0)
      panic("Overlap in ndx file!\n start = %d frames, seen %d frames",
	    ndx_sec->get_start_time(), total_frames_seen());
    if(nskip > 0) {
      INT noskip = dictionary->skip_frames(nskip);
      if(noskip > 0)
	panic("Failed to skip %d frames (skipped %d)\n", nskip, nskip - noskip);
    } 
    frame_offset = ndx_sec->get_start_time();
    set_sentence_length(ndx_sec->get_end_time() - frame_offset);
    if(Debug_flag.ndx){
      printf("NDX[%s]: start_time = %d (%f)  end_time = %d (%f)\n",
	     ndx_sec->get_id(),
	     ndx_sec->get_start_time(), ndx_sec->get_float_start_time(),
	     ndx_sec->get_end_time(), ndx_sec->get_float_end_time());
      printf("Skipping %d frames\n", nskip);
      printf("Frame offset = %d,  sentence_length = %d\n", 
	     frame_offset, sentence_length);
      fflush(stdout);
    }
  } else if(inc_offset) {
    frame_offset = total_frames_seen();
  } else {
    frame_offset = 0;
  }

//  reset_outputs = TRUE;
//  restarted = FALSE;
//  do {  // loop over possible restarts with increased beamwidths
//    restart = FALSE;
  stacks->clear();
  //  new_sentence(reset_outputs, !ndx_output);  // sentence len already set
  new_sentence(TRUE, !ndx_output);  // sentence len already set


  language->reset();

  //  if(reset_outputs) 
  dictionary->init_lub(0);

  if(sentence_length == 0) {
    delete hlist;
    delete delhlist;
    return NULL;
  }

  nhyps_created = 0;
  nhyps_used = 0;
  nnodes_used_sum = 0;
  nnodes_used_sum2 = 0;
  nnodes_pw = 0;
  nlatt_nodes = 0;
  firsth = new HYP(Sentence_start_index, 0, 0, 0, NULL, Sentence_start_index);
  firsth_start = 0;
  assert(firsth->get_total_hyps() == 1);
  pq = new STACK();
  pq->insert(firsth);
  stacks->enqueue(pq);
  
  if(Debug_flag.scheduler)
    printf("Decoding sentence %d (length = %d)\n", sentence_count, sentence_length);

  if(online_output)
    init_online_output();

  if(pipe_lattice){
    LATTICE_LINK::nlattice_elt = 0;
    printf("# Sentence %d\n", sentence_count);
  }

  if(build_lattice) {
    while(lattices->size > 0){
      LIST<LATTICE_LINK*>* ll = lattices->pop();
      while(ll->size > 0)
	delete ll->pop();
      delete ll;
    }	    
    while(silends->size > 0)
      delete silends->pop();
  }

  for(start_time = 0; start_time < sentence_length; start_time++){
    //      fprintf(stderr, "%d total hypotheses\n", firsth->get_total_hyps());
    //    if (online && start_time > 1) {
    //      delete[] outputs->dequeue();
    //      lub->dequeue();
    //      def_lub_inc->dequeue();
    //    }


    if(build_lattice){
      lattices->push(new LIST<LATTICE_LINK*>());
      silends->push(new LIST<INT>());
    }

    assert(hlist->size == 0);
    assert(delhlist->size == 0);

    if(stacks->size == 0){
      //      if(online){
      //	if(Verbose)
      //	  printf("Cannot restart in online mode! Moving to next utterance.\n");
      //	dictionary->eat_acoustic_input(start_time);
      //	return(new STACK());
      //      }
      // Alternatively could increase beam and try again...
      //      assert(smask <= maxmask);
      //      reset_outputs = FALSE;
      //      beam += (INT)log_scale;
      //      state_beam += (INT)log_scale;
      //      restart = TRUE;
      //      restarted = TRUE;
      //      sentence_count |= smask;
      //      smask = smask << 1;
      //      if(Verbose) {
      //	printf("Sentence %d: Ran out of hyps at frame %d\n", 
      //	       sentence_count & sunmask, start_time);
	//	printf("Temporarily increasing beamwidth:  beam = %d  state_beam = %d\n",
	//	       beam, state_beam);
      //	fflush(stdout);
      //      }
      //      break;
      if(smask > maxmask)
	panic("Too many beam increments!");
      sentence_count |= smask;
      smask = smask << 1;
      hlist->push(firsth);
      beam += (INT)log_scale;
      state_beam += (INT)log_scale;
      start_time = firsth->get_time();
      if(Verbose) {
	if(online_output) printf("\n");
	printf("Increasing beam = %.1f  state_beam = %.1f\n", 
	       logdecodelog(beam), logdecodelog(state_beam));
	printf("Rewinding to start_time = %d\n", start_time);
	fflush(stdout);
      }
    } else {
      if(Debug_flag.stack)
	printf("start_time = %d,  %d stacks,  ", start_time, stacks->size);

      // Dequeue the stack for this time
      pq = stacks->dequeue();
      if(Debug_flag.stack) printf("[%d] ", pq->size);
      
      // No hypotheses ending at this time
      if(pq->size == 0){
	if(Debug_flag.stack) printf("0 hypotheses.\n");
	delete pq;
	continue;
      }
      // Make a list of the top nhyp hypotheses
      if(pq->top()->get_time() != start_time)
	panic("t = %d, hypothesis time = %d\n", start_time, h->get_time());



      // Delete hypotheses ending with Sentence_end
      if(!(Sentence_start_index == Sentence_end_index || Cross_sentence)) {
	for(INT i = 0; i < pq->size; i++) {
	  if(pq->get(i)->get_word() == Sentence_end_index)
	    delhlist->push(pq->remove(i));
	}
      }

      // Delete those hyps that are out of beam 
      //  (min-heap so the worst are at the top) 
      if(start_time < proc_frames_seen() && start_time > 0){
	while(pq->size > 0) {
	  h = pq->top();
	  if(lub->get(start_time-1) - h->get_logprob() > beam)
	    delhlist->push(pq->pop());
	  else
	    break;
	}
      }

      garbage_collect(delhlist);

      // If no valid hypotheses
      if(pq->size == 0) {
	// No hypotheses within the beam
	if(Debug_flag.stack) printf("0 hypotheses\n");
	delete pq;
	continue;
      }

      while(pq->size > 0) {
	h = pq->pop();
	hlist->push(h);
      }
      delete pq;

    }

    nnodes_used = 0;
    


    // most probable hypothesis
    top_prob = hlist->top()->get_logprob();
    hlist->init_cursor();
    while(hlist->next())
      hlist->read()->compute_delta(top_prob);

    // incremental output options
    if((Debug_flag.how_are_we_doing && start_time%25==0) ||
       (inc_output > 0 && start_time%inc_output==0)) {
      hlist->top()->display(dictionary->vocabulary, gstdout, start_time);
      fflush(stdout);
    }

    // dump out lattice elements
    if(build_lattice && start_time > 0){
      LIST<WORD_TIME_PAIR*>* wordlist = new LIST<WORD_TIME_PAIR*>();
      LATTICE_LINK *le;
      INT lestart;
      nlatt_nodes++;
      hlist->init_cursor();
      while(hlist->next()) {
	h = hlist->read();
	if(Debug_flag.lattice)
	  h->log_display(dictionary->vocabulary, gstdout);
	// Always have at least <s> before 1st word 
	if(h->get_previous() != NULL && 
	   h->get_previous()->get_word() == Silence_index) {
	  // write interword-pause element
	  HYP* silh = h->get_previous();
	  INT pet = silh->get_time();
	  INT pp = silh->get_word_logprob();
	  // Always have at least <s> before <SIL>
	  lestart = silh->get_previous()->get_time();
	  if(Debug_flag.lattice)
	    printf("silends->get(%d)->push(%d) = %d\n", lestart, pet, silends->get(lestart)->index(pet));
	  if(silends->get(lestart)->index(pet) == -1) {
	    silends->get(lestart)->push(pet);
	    le = new LATTICE_LINK(lestart, pet, Silence_index, pp);
	    if(pipe_lattice) {
	      le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
	      delete le;
	    } else{
	      if(Debug_flag.lattice)
		le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
	      lattices->get(lestart)->push(le);
	    }
	  }
	}
	INT wd = h->get_word();
	lestart = h->get_previous()->get_time();
	WORD_TIME_PAIR* wtp = new WORD_TIME_PAIR(wd, lestart);
	BOOL new_link = TRUE;
	wordlist->init_cursor();
	while(wordlist->next())
	  if(*(wordlist->read()) == *wtp){
	    new_link = FALSE;
	    break;
	  }
	if(new_link) {
	  wordlist->push(wtp);
	  le = new LATTICE_LINK(lestart,
				h->get_time(),
				wd,
				h->get_word_logprob(),
				h->get_pron_version(),
				h->get_pron_prob());
	  if(pipe_lattice) {
	    le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
	    delete le;
	  } else{
	    if(Debug_flag.lattice)
	      le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
	    lattices->get(lestart)->push(le);
	  }
	}
      }
      delete wordlist;
    }

    if(cache_lm && language->get_order() > 1) {
      language->update_cache(hlist);
    }

    if(Debug_flag.stack) {
      if(start_time < proc_frames_seen() && start_time > 0) {
	dl = lub->get(start_time-1) - hlist->top()->get_logprob();
	printf("Extending %d hypotheses\n(lub = %d, delta-lub = %d)\n", 
	       hlist->size, lub->get(start_time-1), dl);
      } else {
	printf("Extending %d hypotheses\n", hlist->size);
      }
      INT i = 0;
      hlist->init_cursor();
      while(hlist->next()){
	hlist->read()->log_display(dictionary->vocabulary, gstdout); 
	if (++i >= MIN(n_hyps_db, hlist->size))
	  break;
      }
      fflush(stdout);
    }

    if(Debug_flag.decode_stats)
      Decode_inf_file->fprintf("%d ", hlist->size);

    nhyps_used += hlist->size;

    // This extends the current set of hypotheses, writing new ones into
    //  the appropriate stack
    if(Aligning)
      nh = dictionary->align(start_time, hlist, language, stacks);
    else
      nh = dictionary->extend(start_time, hlist, language, stacks);

    if(start_time == 0 || hlist->top() == firsth) 
      hlist->pop(); // don't delete firsth
    else
      garbage_collect(hlist);

    garbage_collect(Hyp_garbage_list);

    if(online_output)
      hyp_chain_forward();
	
    if(Debug_flag.stack) {
      printf("Created %d new hypotheses\n", nh);
      fflush(stdout);
    }
    nhyps_created += nh;
    max_nhyps_created = MAX(max_nhyps_created, nh);
    nnodes_used_sum += dictionary->node_count;
    nnodes_used_sum2 += nnodes_used;
    max_nnodes_used = MAX(max_nnodes_used, dictionary->node_count);
    max_nnodes_used2 = MAX(max_nnodes_used2, nnodes_used);
  }

  //  } while(restart);
  //  if(restarted) {
  if(beam > save_beam) {
    sentence_count &= sunmask;
    beam = save_beam;
    state_beam = save_sbeam;
    if(Verbose) {
      printf("Resetting beam = %.1f  state_beam = %.1f\n", 
	     logdecodelog(beam), logdecodelog(state_beam));
      fflush(stdout);
    }
  }
  //  }

  if(sentence_length == 0){
    // end of all sentences
    res = NULL;
  }  else {
    res = new STACK();
    if (stacks->size == 0) {
      fprintf(stderr, "No final stack!\n");
    } else {
      assert(stacks->size == 1);
      pq = stacks->dequeue();
      if(pq->size == 0){
	fprintf(stderr, "No final hypothesis!\n");
      } else {
	if(Debug_flag.stack)
	  printf(" %d final hypotheses\n", pq->size);
	// dump out lattice elements from final stack
	if(build_lattice){
	  LIST<WORD_TIME_PAIR*>* wordlist = new LIST<WORD_TIME_PAIR*>();
	  LATTICE_LINK *le;
	  INT lestart;
	  nlatt_nodes++;
	  while(pq->size > 0){
	    h = pq->pop();
	    if(Debug_flag.lattice)
	      h->log_display(dictionary->vocabulary, gstdout);
	    // Always have at least <s> before 1st word 
	    HYP* silh = h->get_previous();
	    if(silh != NULL && silh->get_word() == Silence_index) {
	      // write interword-pause element
	      INT pet = silh->get_time();
	      INT pp = silh->get_word_logprob();
	      // Always have at least <s> before <SIL>
	      lestart = silh->get_previous()->get_time();
	      if(silends->get(lestart)->index(pet) == -1) {
		silends->get(lestart)->push(pet);
		le = new LATTICE_LINK(lestart, pet, Silence_index, pp);
		if(pipe_lattice) {
		  le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
		  delete le;
		} else{
		  if(Debug_flag.lattice)
		    le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
		  lattices->get(lestart)->push(le);
		}
	      }
	    }
	    INT wd = h->get_word();
	    if(wd == Sentence_end_index)
	      wd = Silence_index;
	    lestart = h->get_previous()->get_time();
	    WORD_TIME_PAIR* wtp = new WORD_TIME_PAIR(wd, lestart);
	    if(lestart < h->get_time()) {
	      BOOL new_link = TRUE;
	      wordlist->init_cursor();
	      while(wordlist->next())
		if(*(wordlist->read()) == *wtp){
		  new_link = FALSE;
		  break;
		}
	      if(new_link) {
		wordlist->push(wtp);
		le = new LATTICE_LINK(lestart,
				      h->get_time(),
				      wd,
				      h->get_word_logprob(),
				      h->get_pron_version(),
				      h->get_pron_prob());
		if(pipe_lattice) { 
		  le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
		  delete le;
		} else{
		  if(Debug_flag.lattice)
		    le->print_with_frames(gstdout, dictionary->vocabulary, le->num_links()-1);
		  lattices->get(lestart)->push(le);
		}
	      }
	    }
	    if(h->get_word() != Sentence_end_index){
	      HYP* xh = new HYP(Sentence_end_index, h->get_time(),
			h->get_logprob() + language->end_of_sentence_prob(h), 0, h);
	      res->insert(xh);
	    } else {
	      res->insert(h);
	    }
	  }
	  delete wordlist;
	} else {
	  HYP *hyp = NULL;
	  while(pq->size > 0){
	    hyp= pq->pop();
	    if(hyp->get_word() != Sentence_end_index){
	      HYP* xh = new HYP(Sentence_end_index, hyp->get_time(),
		   hyp->get_logprob() + language->end_of_sentence_prob(hyp), 0, hyp);
	      res->insert(xh);
	    } else {
	      res->insert(hyp);
	    }
	  }
	}
      }
      delete pq;
    }


    if(Debug_flag.stack)
      printf(" %d final hypotheses\n", res->size);

    if(online_output) {
      while(res->size > 0)
	hlist->push(res->pop());

      if(hlist->size > 0) {
	HYP *toph = hlist->pop();
	garbage_collect(hlist);
	// the only HYPs left should be the ones that need printing
	hyp_chain_forward();
	// now only toph should be left - deal with this explicitly
	if(firsth != toph || firsth->get_total_hyps() != 1){
	  firsth->info(dictionary->vocabulary,dictionary->phonenames);
	  toph->info(dictionary->vocabulary,dictionary->phonenames);
	  panic("firsth = %p toph = %p,  total_hyps = %d\n",
		firsth, toph, firsth->get_total_hyps());
	}

	// print out the final word
	firsth->online_display(dictionary->vocabulary, output_fp);
	if(state_decode_output) {
	  firsth->online_state_display(dictionary->phonenames, firsth_start, 
				       frame_offset, frame_shift, state_dec_fp);
	}
	if(ctm_output){
	  if(Debug_flag.ndx) 
	    ctm_fp->fprintf("**%.3f** ", 
			    (FLOAT)(start_time-firsth->get_time())*frame_shift);
	  firsth->online_ctm_display(dictionary->vocabulary, firsth_start, 
				     frame_offset, frame_shift, ctm_fp, 
				     sentence_count, 'A');
	}
	if(srt_output)
	  firsth->online_srt_display(dictionary->vocabulary, firsth_start, 
				     frame_offset, frame_shift, srt_fp);
	if(Verbose) 
	  firsth->online_log_display(dictionary->vocabulary, gstdout);	
	gstdout->fprintf("\n");
      }
    }

    if(Verbose){
      printf("Created %d hypotheses, extended %d\n", nhyps_created, nhyps_used);
      printf("Activated %d nodes (%.1f per start)\n", nnodes_used_sum2,
	     (FLOAT)nnodes_used_sum2/(FLOAT)sentence_length);
      printf("Propagated %d times: %.1f nodes per frame  (%.1f frames per node)\n",
	     nnodes_used_sum, (FLOAT)nnodes_used_sum/(FLOAT)sentence_length,
	     (FLOAT)nnodes_used_sum/(FLOAT)nnodes_used_sum2); 
      if(Debug_flag.tree_stats)
	printf(" (%d above pw threshold)\n", nnodes_pw);
      //      printf("Reactivated %d nodes out of order (%.1f per frame)\n",
      //	     nreactivated, (FLOAT)nreactivated/(FLOAT)sentence_length);
      printf("Pruned %d phones (%.1f%%) by PDP\n", nphones_pruned, 
	     (100.0 * nphones_pruned) / (FLOAT)(sentence_length*nout));
#ifdef LM_STATS
      language->access_info();
#endif
      fflush(stdout);
      total_hyps_created += nhyps_created;
      total_hyps_used += nhyps_used;
      total_nodes_used += nnodes_used_sum;
      total_nodes_used2 += nnodes_used_sum2;
      total_phones_pruned += nphones_pruned;
      if(Debug_flag.tree_stats)
	total_nodes_pw += nnodes_pw;
    }
  }
  delete firsth;
  delete hlist;
  delete delhlist;
  return res;
}

VOID
SCHEDULER::build(){
  GENFILE *phone_fp = NULL;
  GENFILE *prior_fp = NULL;
  GENFILE *lex_fp = NULL;
  GENFILE* lm_fp = NULL;
  char *lm_fname = NULL;
  INT nwords = 0;

  // Make the lookup table for logadd
  DOUBLE diff;
  for(INT j = 0; j < logadd_lookup_length; j++) {
    diff = (DOUBLE)j * one_over_log_scale;
    logadd_lookup[j] = (INT)(log10(1.0 + exp10(-diff))*log_scale + 0.5);
  }

  if(param->get_bool("output_file")){
    output_fp = new GENFILE(param->get_string("output_file", "OP"), "w");
  }  else {
    output_fp = gstdout;
  }

  inc_output = 0;
  //  if(param->get_bool("demo"))
  //     inc_output = 1;
  if(param->get_bool("inc_output"))
    inc_output = param->get_int("inc_output", 25);
  
  Sentence_start_symbol = strdup(param->get_string("sentence_start_symbol", "<s>"));
  Sentence_end_symbol = strdup(param->get_string("sentence_end_symbol", "</s>"));
  Paragraph_symbol = strdup(param->get_string("paragraph_symbol", "<p>"));
  Article_symbol = strdup(param->get_string("article_symbol", "<art>"));
  Silence_symbol = strdup(param->get_string("silence_symbol", "<SIL>"));
  Unknown_word_symbol = strdup(param->get_string("unknown_word_symbol", "<UNK>"));
  Cross_sentence = param->get_bool("cross_sentence");

  State_decode = param->get_bool("phone_decode");
  if(State_decode) {
    char *sdnm = strdup(param->get_string("phone_decode", "phone.decode"));
    state_dec_fp = new GENFILE(sdnm, "w");
    delete sdnm;
    state_decode_output = TRUE;
  }  

  dump_lub = param->get_bool("dump_lub");
  if(dump_lub){
    dump_lub_name = strdup(param->get_string("dump_lub", "LUB"));
    State_decode = TRUE;
  }

  build_lattice = param->get_bool("lattice");
  if(build_lattice){
    lattice_name = strdup(param->get_string("lattice", "LATTICE"));
  }

  //  build_detailed_lattice = param->get_bool("detailed_lattice");
  //  if(build_detailed_lattice){
  //    detailed_lattice_name = strdup(param->get_string("detailed_lattice", "DETAILED_LATTICE"));
  //  }

  // DCA 18/MAR/96: pipe lattice option
  pipe_lattice = param->get_bool("pipe_lattice");
  if(pipe_lattice) {
    if(build_lattice)
      fprintf(stderr, "Warning: piping lattice not writing to file\n(both -lattice and -pipe_lattice specified)\n");
    build_lattice = TRUE;
    LATTICE_LINK::pipe = TRUE;
  }

  if(build_lattice)
    Pause_decode = TRUE;
  
  //  if(pipe_lattice){ 
  //     if((build_lattice && build_detailed_lattice) || !(build_lattice || build_detailed_lattice)){
  //        panic("-pipe_lattice can only be used in conjunction with EITHER -build_lattice OR -build_detailed_lattice\n");
  //     }

  // DCA 25/MAR/96: ctm output option
  frame_shift = 0.001 * (FLOAT)param->get_int("frame_shift", 16);

  ndx_output = param->get_bool("ndx");
  if(ndx_output){
    char *ndx_nm = strdup(param->get_string("ndx", "ndx"));
    GENFILE *ndx_fp = new GENFILE(ndx_nm, "r");
    delete ndx_nm;
    Pause_decode = TRUE;
    init_ndx(ndx_fp);
    delete ndx_fp;
  }

  ctm_output = param->get_bool("ctm");
  if(ctm_output){
    ctm_filename = strdup(param->get_string("ctm", "CTM", 0));
    Pause_decode = TRUE;
    if(ndx_output){
      ctm_fp = new GENFILE(ctm_filename, "w");
      ctm_fp->fprintf("; %s\n", ndx_hdr);
    }      
  }

  srt_output = param->get_bool("srt");
  if(srt_output){
    srt_filename = strdup(param->get_string("srt", "SRT", 0));
    Pause_decode = TRUE;
    if(ndx_output){
      srt_fp = new GENFILE(srt_filename, "w");
      srt_fp->fprintf("%s\n", ndx_hdr);
    }      
  }

  online = param->get_bool("online");
  if(online && !(build_lattice || dump_lub))
    online_output = TRUE;

  if((ctm_output || srt_output || online) && !ndx_output &&
     !param->get_bool("reset_time_each_utt"))
    inc_offset = TRUE;
  else
    inc_offset = FALSE;


  Forward_process = param->get_bool("forward_process");
  if(Forward_process && State_decode) {
    fprintf(stderr, "-phone_decode not compatible with -forward_process, using Viterbi\n");
    Forward_process = FALSE;
  }

  if(param->get_bool("new_lub") || param->get_bool("demo")) {
    strict_lub = TRUE;
    if(Verbose && param->get_bool("demo"))
      printf("-demo option setting -new_lub\n");
  }  else {
    strict_lub = FALSE;
  }
  //  if(strict_lub && Forward_process) {
  //    fprintf(stderr, "-new_lub and -forward_process not compatible\n");
  //    fprintf(stderr, "Ignoring -forward_process");
  //    Forward_process = FALSE;
  //  }

  if(param->get_bool("align")){
     Aligning = TRUE;
     State_decode = TRUE;
  }

  if(strict_lub || State_decode || Pause_decode)
    Backtrace = TRUE;


  // statics
  pathword_threshold = param->get_int("pathword_threshold", 0);
  FLOAT beam_default, state_beam_default;
  // demo
  beam_default = 4.; 
  state_beam_default = 5.;
  if(Verbose && param->get_bool("demo"))
    printf("-demo option setting -beam %g -state_beam %g\n", 
	   beam_default, state_beam_default);
  beam = (INT)(param->get_float("beam", beam_default) * (FLOAT)log_scale);
  state_beam = (INT)(param->get_float("state_beam", state_beam_default) *
		     (FLOAT)log_scale);     

  if(Aligning)
    dictionary = new LINEAR_LEXICON();
  else
    dictionary = new TREE_LEXICON();

  if(Debug_flag.scheduler) {
    fprintf(stderr, "Created dictionary\n");
  }
  
  cache_lm = TRUE;
  if (param->get_bool("tagprobs")) {
    lm_fname = strdup(param->get_string("ngram", "stdin"));
    language = new TAGGED_NGRAM();
    cache_lm = FALSE;
    // hack! so we don't need use_unk
    delete (char*)Unknown_word_symbol;
    Unknown_word_symbol = strdup("<REAL_UNKNOWN>");
  } else if(param->get_bool("blind_mixture_ngram")) {
    lm_fname = strdup(param->get_string("blind_mixture_ngram", "stdin"));
    language = new BLIND_MIXTURE_NGRAM();
    cache_lm = FALSE;
  } else if(param->get_bool("lsa_mixture_ngram")) {
    lm_fname = strdup(param->get_string("lsa_mixture_ngram", "stdin"));
    language = new LSA_MIXTURE_NGRAM();
    cache_lm = FALSE;
  } else if(param->get_bool("ngram")) {     
    lm_fname = strdup(param->get_string("ngram", "stdin"));
    language = new NGRAM();
  } else if(param->get_bool("bigram")){
    lm_fname = strdup(param->get_string("bigram", "stdin"));
    language = new BIGRAM();
  } else if(param->get_bool("trigram")){
    lm_fname = strdup(param->get_string("trigram", "stdin"));
    language = new TRIGRAM();
  } else if(param->get_bool("nogram")) {
    lm_fname = NULL;
    language = new NOGRAM();
    cache_lm = FALSE;
  } else if(Aligning){
    language = new ALIGN_LM();
    cache_lm = FALSE;
    lm_fname = strdup(param->get_string("align", "REFS"));
  } else {
    panic("Need to specify -ngram, -trigram, -bigram -align or -nogram");
  }

  if(param->get_bool("no_lm_cache"))
    cache_lm = FALSE;

  if(Debug_flag.scheduler) {
    fprintf(stderr, "Created LM\n");
  }
  
  // Use probs of <UNK> in LM?
  language->use_unk = (param->get_bool("use_unk") || param->get_bool("write_lm"))
                      && !language->tagged;

  FLOAT uf = param->get_float("unk_factor", 1.0);
  if(uf < 0.0 || uf > 1.0)
    panic("0.0 < unk_factor <= 1.0 !");
  if(uf == 0.0){
    language->unk_factor = 0xffff;
    language->use_unk = FALSE;
  }
  else if(uf == 1.0)
    language->unk_factor = 0;
  else
    language->unk_factor = (USHORT)MIN(50000,-logcode(uf));
  // Don't auto set to 1/nunigrams
    //    language->unk_factor = 0xffff;

  average_lm_contexts = param->get_int("average_lm_contexts", 0);

  lmscale = param->get_int("lm_scale", 1);
  lmoffset = param->get_int("lm_offset", 0);

  if(param->get_bool("no_decode") || param->get_bool("text_decode") ||
      !(param->get_bool("phone_models") || param->get_bool("dictionary")))
    cache_lm = FALSE;

  if(Verbose){
    if(cache_lm)
      printf("LM Caching\n");
    else
      printf("LM Caching Off\n");
  }

  if(Debug_flag.scheduler) {
    fprintf(stderr, "Ready to build language\n");
  }

  if(!strcmp(lm_fname, "stdin")){
    if(param->get_bool("no_decode"))
      lm_fp = gstdin;
    else
      panic("Can only have LM file piped in using '-no_decode'");
  } else {
    if(lm_fname != NULL)
      lm_fp = new GENFILE(lm_fname,"r");
    else
      lm_fp = NULL;
  }

  // set up the grammar
  // this also writes the dictionary->{vocabulary,vocab_hash}
  nwords = language->build(lm_fp, dictionary, lmscale, lmoffset, cache_lm,
			     maxhyps);

  delete lm_fp;

  HYP::mixture_lm_size = language->get_ncomponents();
  if(Debug_flag.mixlm)
    printf("HYP::mixture_lm_size = %d\n", HYP::mixture_lm_size);
  HYP::lm_context_length = language->get_order();
  HYP::output_pron_version = param->get_bool("output_pron_version");
  HYP::output_silence = param->get_bool("output_silence");


  if(language->tagged){
    GENFILE* tag_fp = new GENFILE(param->get_string("tagprobs", "TAGS"), "r");
    nwords = ((TAGGED_NGRAM*)language)->build_tags(tag_fp, dictionary);
    delete tag_fp;
  }


  if(language->use_unk == TRUE && Unknown_word_index == -1) {
    fprintf(stderr, "WARNING:  <UNK> is not in the LM, but -use_unk specified!\n");
    fprintf(stderr, " (Expect assert error if pronunciation tree contains words not in LM)\n");
  }

  if(Cross_sentence && Sentence_start_index != Sentence_end_index) {
    Cross_sentence = FALSE;
    gstderr->fprintf("No cross sentence decoding when Sentence_start_symbol != Sentence_end_symbol\n");
  }

  if(param->get_bool("write_lm")) {
    STRING wlm_fname = param->get_string("write_lm", "lm.bin");
    lm_fp = new GENFILE(wlm_fname, "w");
    language->write(lm_fp, dictionary);
    delete lm_fp;
    if(Verbose)
      printf("Wrote binary LM file %s\n", wlm_fname);
  }
  
  if(param->get_bool("print_unigrams")) {
    GENFILE *fp = NULL;
    if(param->get_size("print_unigrams") == 0)
      fp = gstdout;
    else
      fp = new GENFILE(param->get_string("print_unigrams"), "w");
    language->print_unigrams(fp, dictionary);
    delete fp;
  }

  if(param->get_bool("print_lm")) {
    GENFILE *fp = NULL;
    if(param->get_size("print_lm") == 0)
      fp = gstdout;
    else
      fp = new GENFILE(param->get_string("print_lm"), "w");
    language->print(fp, dictionary);
    delete fp;
  }

  if(param->get_bool("text_decode")) {
    GENFILE *text_fp;
    BOOL ignore_nl = TRUE;
    if(param->get_bool("nl_is_end_sent"))
      ignore_nl = FALSE;
    STRING text_fname = param->get_string("text_decode", "stdin");
    if(!strcmp(text_fname, "stdin")) {
      if(!strcmp(lm_fname, "stdin"))
	panic("Cannot have both LM and text input on stdin");
      text_fp = gstdin;
    } else {
      text_fp = new GENFILE(text_fname, "r");
    }
    DOUBLE noUNKperp = 0.0;
    DOUBLE perp = language->text_decode(text_fp, output_fp, dictionary, noUNKperp, ignore_nl);    
    if(language->use_unk)
      printf("Test Set Perplexity = %f (%f including <UNK>)\n", noUNKperp, perp);
    else
      printf("Test Set Perplexity = %f\n", perp);
    fflush(stdout);
  }
  delete lm_fname;

  // exit if no decoding is required
  if(param->get_bool("no_decode") || param->get_bool("text_decode") ||
      !(param->get_bool("phone_models") || param->get_bool("dictionary")))
    exit(0);
  

  // Read in phone set and pronunciation dictionary

  // Open phone and prior files
  // Either 1.  -phi <file>
  //   Or   2.  -phone_models <file1> -priors <file2> 
  BOOL phiflag = FALSE;
  if(param->get_bool("phone_models") || param->get_bool("phi")) {
    phiflag = param->get_bool("phi");
    if(param->get_bool("phone_models") && phiflag) 
      panic("Only specify one of -phone_models or -phi\n");

    if(phiflag) {
      if(param->get_bool("priors"))
	fprintf(stderr, "-priors ignored when -phi specified");
      phone_fp = new GENFILE(param->get_string("phi", "PHI"), "r");
      prior_fp = NULL;
    } else {
      if(!param->get_bool("priors"))
      fprintf(stderr, "must specify -priors with -phone_models");
      phone_fp = new GENFILE(param->get_string("phone_models", "PHONE"), "r");
      prior_fp = new GENFILE(param->get_string("priors", "PRIORS"), "r");
    }
  } else {
    panic("Must specify '-phone_models' and '-prior' OR '-phi'\n");
  }

  // Open dictionary file
  if(!param->get_bool("dictionary"))
    panic("Must specify '-dictionary <file>'\n");
  lex_fp = new GENFILE(param->get_string("dictionary", "DICTIONARY"), "r");

  // Builds dictionary and phone set and initialises acoustic input
  //  if -no_decode is not set
  dictionary->build(param, nwords, phiflag, lex_fp, phone_fp, prior_fp);

  delete phone_fp;
  delete prior_fp;
  delete lex_fp;

  if(param->get_bool("no_decode"))
    exit(0);


  //  dictionary->set_acoustic_input(param);

  stacks = new QUEUE<STACK*>(500);
  Hyp_garbage_list = new LIST<HYP*>(128);

  //  if(Debug_flag.tree_stats)
  //    for(INT i = 0; i < dictionary->depth_count->size; i++)
  //      depth_count->enqueue(0);

  return;
}

VOID
SCHEDULER::parse(INT argc, STRING argv[]){
  STRING param_files = "parameter_files";

  if (argc == 1) {
    param->help_message();
    exit(0);
  }

  // The parameter list consists of a sequence of parameter names
  // and values.  Parameter names always start with a minus. 
  // Values can be numbers or strings (not starting with minus).
  // Any number of values (including 0) can follow a parameter name.
  // Access to these values is through functions such as param->get_int("name"),  
  // param->get_float("name",default) and param->get_string("name",default,index)

  // Set the default parameter name for values that occur before any
  // parameter name to "parameter_files".
  param->set_default_param_name(param_files);

  // Add the parameters from the command line to the parameter list.
  param->parse_cmd_line(argc,argv,OVERRIDE);

  // Any values that come before the first parameter name are added as
  // "parameter_files".  In this case, we use these arguments 
  // as file names that contain more parameters.
  for(INT i = 0; param->get_string(param_files,NULL,i) != NULL; i++) {
    param->parse_file(param->get_string(param_files,NULL,i),UNDERRIDE);
  }

  if(param->get_bool("version"))
    exit(0);

  if(param->get_bool("help")) {
    param->help_message();
    exit(0);
  }

  if (param->get_bool("d")) {
    int i;
    if (param->get_size("d") == 0)
      Debug_flag.set_all();
    else {
      for(i = 0; i < param->get_size("d"); i++) {
        Debug_flag.set((CCHAR*)param->get_string("d","?",i));
      }
    }
  }

  Verbose = param->get_bool("verbose");
}

VOID
SCHEDULER::init_online_output(){
  if(ndx_output) {
    if(srt_output) ndx_sec->print_start(srt_fp);
    if(ctm_output) ndx_sec->print_start(ctm_fp, ";");
  } else {
    if(srt_output) {
      CHAR srtnm[1024];
      CHAR strbuf[BUFSIZ];
      if ((strcmp(&srt_filename[strlen(srt_filename) - 3], ".gz") == 0) ) {
	strncpy(strbuf, srt_filename, strlen(srt_filename) - 3);
	strbuf[strlen(srt_filename)-3] = '\0';
	sprintf(srtnm, "%s.%d.gz", strbuf, sentence_count);
      } else if ((strcmp(&srt_filename[strlen(srt_filename) - 2], ".Z") == 0) ) {
	strncpy(strbuf, srt_filename, strlen(srt_filename) - 2);
	strbuf[strlen(srt_filename)-2] = '\0';
	sprintf(srtnm, "%s.%d.Z", strbuf, sentence_count);
      } else {
	sprintf(srtnm, "%s.%d", srt_filename, sentence_count);
      }
      srt_fp = new GENFILE(srtnm, "w");
    }
    if(ctm_output) {
      CHAR ctmnm[1024];
      CHAR strbuf[BUFSIZ];
      if ((strcmp(&ctm_filename[strlen(ctm_filename) - 3], ".gz") == 0) ) {
	strncpy(strbuf, ctm_filename, strlen(ctm_filename) - 3);
	strbuf[strlen(ctm_filename)-3] = '\0';
	sprintf(ctmnm, "%s.%d.gz", strbuf, sentence_count);
      } else if ((strcmp(&ctm_filename[strlen(ctm_filename) - 2], ".Z") == 0) ) {
	strncpy(strbuf, ctm_filename, strlen(ctm_filename) - 2);
	strbuf[strlen(ctm_filename)-2] = '\0';
	sprintf(ctmnm, "%s.%d.Z", strbuf, sentence_count);
      } else {
	sprintf(ctmnm, "%s.%d", ctm_filename, sentence_count);
      }
      ctm_fp = new GENFILE(ctmnm, "w");
    }
  }
}


#define LATBUFSIZ 125
VOID
SCHEDULER::output_hyps(STACK *hyps){
  CHAR latnm[1024];
  CHAR dumpnm[1024];
  CHAR ctmnm[1024];
  CHAR srtnm[1024];
  LIST<HYP*>* hlist = new LIST<HYP*>();

  if(!online_output) {
    while(hyps->size > 0)
      hlist->push(hyps->pop());
  
    if(Verbose){
      printf("%d final hypotheses\n", hlist->size);
      fflush(stdout);
    }

    if(hlist->size == 0){
      delete hlist;
      return;
    }
  }

  if(ndx_output && sentence_length != ndx_sec->get_end_time() - frame_offset){
    fprintf(stderr, "%d frames in ndx section (expecting %d)\n",
	    sentence_length, ndx_sec->get_end_time() - frame_offset);
    ndx_sec->info(gstderr);
  }


  if(online_output) {
    // deal with the ouput files at the utterance end
    if(ndx_output) {
      if(ctm_output) ndx_sec->print_end(ctm_fp, ";");
      if(srt_output) ndx_sec->print_end(srt_fp);
    } else {
      if(ctm_output) delete ctm_fp;
      if(srt_output) delete srt_fp;
    }
    output_fp->fprintf("\n");
    output_fp->fflush();
    if(Verbose)
      printf("\n");    // log file

  } else if(hlist->size == 0){
    output_fp->fprintf("%d <NULL>\n", sentence_count);

  } else {
    HYP *toph = hlist->top();

    toph->display(dictionary->vocabulary, output_fp, sentence_count,
		  output_logprob);
    output_fp->fflush();

    if(state_decode_output) {
      toph->display(dictionary->vocabulary, state_dec_fp, sentence_count,
		    output_logprob);
      toph->state_display(dictionary->phonenames, frame_offset, frame_shift, 
			  state_dec_fp);
    }    

    if(ctm_output) {
      if(ndx_output) {
	ndx_sec->print_start(ctm_fp, ";");
	toph->ctm_display(dictionary->vocabulary, frame_offset, frame_shift, ctm_fp, 
			  sentence_count, 'A');
	ndx_sec->print_end(ctm_fp, ";");
	ctm_fp->fflush();
      } else {
	CHAR strbuf[BUFSIZ];
	if ((strcmp(&ctm_filename[strlen(ctm_filename) - 3], ".gz") == 0) ) {
	  strncpy(strbuf, ctm_filename, strlen(ctm_filename) - 3);
	  strbuf[strlen(ctm_filename)-3] = '\0';
	  sprintf(ctmnm, "%s.%d.gz", strbuf, sentence_count);
	} else if ((strcmp(&ctm_filename[strlen(ctm_filename) - 2], ".Z") == 0) ) {
	  strncpy(strbuf, ctm_filename, strlen(ctm_filename) - 2);
	  strbuf[strlen(ctm_filename)-2] = '\0';
	  sprintf(ctmnm, "%s.%d.Z", strbuf, sentence_count);
	} else {
	  sprintf(ctmnm, "%s.%d", ctm_filename, sentence_count);
	}
	ctm_fp = new GENFILE(ctmnm, "w");
	toph->ctm_display(dictionary->vocabulary, frame_offset, frame_shift, ctm_fp, 
			  sentence_count, 'A');
	delete ctm_fp;
      }
    }


    if(srt_output) {
      if(ndx_output) {
	ndx_sec->print_start(srt_fp);
	toph->srt_display(dictionary->vocabulary, frame_offset, frame_shift, srt_fp);
	ndx_sec->print_end(srt_fp);
	srt_fp->fflush();
      } else {
	CHAR strbuf[BUFSIZ];
	if ((strcmp(&srt_filename[strlen(srt_filename) - 3], ".gz") == 0) ) {
	  strncpy(strbuf, srt_filename, strlen(srt_filename) - 3);
	  strbuf[strlen(srt_filename)-3] = '\0';
	  sprintf(srtnm, "%s.%d.gz", strbuf, sentence_count);
	} else if ((strcmp(&srt_filename[strlen(srt_filename) - 2], ".Z") == 0) ) {
	  strncpy(strbuf, srt_filename, strlen(srt_filename) - 2);
	  strbuf[strlen(srt_filename)-2] = '\0';
	  sprintf(srtnm, "%s.%d.Z", strbuf, sentence_count);
	} else {
	  sprintf(srtnm, "%s.%d", srt_filename, sentence_count);
	}
	srt_fp = new GENFILE(srtnm, "w");
	toph->srt_display(dictionary->vocabulary, frame_offset, frame_shift, srt_fp);
	delete srt_fp;
      }
    }

    if(dump_lub) {
      LIST<INT> *prob_seq = new LIST<INT>(1024);
      LIST<INT> *pseq2 = new LIST<INT>(1024);
      dictionary->score(toph, language, prob_seq, pseq2);
      assert(prob_seq->size == sentence_length);
      CHAR strbuf[BUFSIZ];
      if ((strcmp(&dump_lub_name[strlen(dump_lub_name) - 3], ".gz") == 0) ) {
	strncpy(strbuf, dump_lub_name, strlen(dump_lub_name) - 3);
	strbuf[strlen(dump_lub_name)-3] = '\0';
	sprintf(dumpnm, "%s.%d.gz", strbuf, sentence_count);
      } else if ((strcmp(&dump_lub_name[strlen(dump_lub_name) - 2], ".Z") == 0) ) {
	strncpy(strbuf, dump_lub_name, strlen(dump_lub_name) - 2);
	strbuf[strlen(dump_lub_name)-2] = '\0';
	sprintf(dumpnm, "%s.%d.Z", strbuf, sentence_count);
      } else {
	sprintf(dumpnm, "%s.%d", dump_lub_name, sentence_count);
      }
    
      GENFILE* dump_fp = new GENFILE(dumpnm, "w");
      
      for(INT i = 0; i < sentence_length; i++){
	dump_fp->fprintf("%d %f %f %f", i, logdecodelog(lub->get(i)), 
			 logdecodelog(prob_seq->get(i)),
			 logdecodelog(pseq2->get(i)));
	//	if(i >= toph->node_time_history()->get(j)-1){
	//	  fprintf(dump_fp, "(%d: %f)", toph->node_time_history()->get(j)-1,
	//		  logdecodelog(toph->node_prob_history()->get(j)));
	//	  j++;
	// 	}
	dump_fp->fprintf("\n");
      }
      printf("\n");
      delete dump_fp;
      delete prob_seq;
      delete pseq2;
    }

    if(Verbose)
      toph->log_display(dictionary->vocabulary, gstdout);	


    if(build_lattice && !pipe_lattice){
      CHAR strbuf[BUFSIZ];
      if ((strcmp(&lattice_name[strlen(lattice_name) - 3], ".gz") == 0) ) {
	strncpy(strbuf, lattice_name, strlen(lattice_name) - 3);
	strbuf[strlen(lattice_name)-3] = '\0';
	sprintf(latnm, "%s.%d.gz", strbuf, sentence_count);
      } else if ((strcmp(&lattice_name[strlen(lattice_name) - 2], ".Z") == 0) ) {
	strncpy(strbuf, lattice_name, strlen(lattice_name) - 2);
	strbuf[strlen(lattice_name)-2] = '\0';
	sprintf(latnm, "%s.%d.Z", strbuf, sentence_count);
      } else {
	sprintf(latnm, "%s.%d", lattice_name, sentence_count);
      }

      GENFILE* lattice_fp = new GENFILE(latnm, "w");

      // clean lattice and build nodes
      LIST<LIST<LATTICE_LINK*>*>* clean_lattice = new LIST<LIST<LATTICE_LINK*>*>();
      LIST<LATTICE_NODE*>* lnl = new LIST<LATTICE_NODE*>();
      LIST<INT>* nodetimes = new LIST<INT>();
      INT nlinks = 0;
      nodetimes->put_and_grow(sentence_length, 0, -1);
      lnl->push(new LATTICE_NODE(sentence_length));

      // backward pass through lattice links
      BOOL new_node_flag = TRUE;
      while(lattices->size > 0){
	LIST<LATTICE_LINK*>* ll = lattices->pop();
	if(ll == NULL) continue;
	LATTICE_LINK *lnk;
	new_node_flag = TRUE;
	LIST<LATTICE_LINK*>* clean_ll = new LIST<LATTICE_LINK*>();
	while(ll->size > 0) {
	  lnk = ll->pop();
	  INT id = nodetimes->get(lnk->get_end_time());
	  if(id == -1){
	    delete lnk;
	    continue;
	  }
	  lnk->put_end_node(id);
	  if(new_node_flag){
	    INT t = lnk->get_start_time();
	    lnl->push(new LATTICE_NODE(t));
	    nodetimes->put(t, lnl->size-1);
	    new_node_flag = FALSE;
	  }
	  lnk->put_start_node(lnl->size-1);
	  clean_ll->push(lnk);
	  nlinks++;
	}
	delete ll;
	if(clean_ll->size > 0)
	  clean_lattice->push(clean_ll);
	else
	  delete clean_ll;
      }
      delete nodetimes;

      // clean up silends
      while(silends->size > 0)
	delete silends->pop();

      assert(nlinks == LATTICE_LINK::nlattice_elt);

      print_lattice_header(lattice_fp, param, toph, dictionary->vocabulary, 
			   sentence_count, nlinks, lnl->size);

      lattice_fp->fprintf("#\n# Node definitions\n#\n");
      INT node_count = 0;
      while(lnl->size > 0) {
	LATTICE_NODE* ln = lnl->pop();
	ln->print(lattice_fp, frame_offset, frame_shift, node_count++);
	delete ln;
      }
      delete lnl;

      lattice_fp->fprintf("#\n# Link definitions\n#\n");
      INT count = 0;
      while(clean_lattice->size > 0) {
	LIST<LATTICE_LINK*> *lat = clean_lattice->pop();
	while(lat->size > 0) {
	  LATTICE_LINK *lnk = lat->pop();
	  lnk->print(lattice_fp, dictionary->vocabulary, count++, node_count);
	  delete lnk;
	}
	delete lat;
      }
      delete clean_lattice;
      delete lattice_fp;
    }
    garbage_collect(hlist);
  }
  gstdout->fflush(); 
  delete hlist;
  delete ndx_sec;
}


// DCA 11/MAR/96: print lattice header
VOID 
SCHEDULER::print_lattice_header(GENFILE *fp, PARAM_TABLE *ptable, 
				HYP *hyp, 
				LIST<STRING> *vocab, INT sc, 
				INT num_links, INT num_nodes){
  fp->fprintf("# Lattice generated by noway version %d.%d", 
	      Version, Subversion);
  if(Patchlevel > 0)
    fp->fprintf("p%d",Patchlevel);
#ifdef TIMING
  fp->fprintf(" %s", ctime(&tt)); // ctime returns a string ending in \n
#else
  fp->fprintf("\n");
#endif
  fp->fprintf("base=%g\n", logbase);
  fp->fprintf("lmscale=%d\n", lmscale);
  fp->fprintf("wdpenalty=%d\n", lmoffset);
  fp->fprintf("amscale=%g\n", amscale);
  ptable->print(fp->get_fp(), "# ");
  fp->fprintf("# Frame shift: %g ms\n", frame_shift*1000.0);
  fp->fprintf("#\n");
  fp->fprintf("# Best hypothesis:\n");
  hyp->display(vocab, fp, sc, TRUE, "# ");
  fp->fprintf("#\n");
  hyp->ctm_display(vocab, frame_offset, frame_shift, fp, 0, '#');
  fp->fprintf("#\n");
  fp->fprintf("# Size line\n");
  fp->fprintf("#\n");
  fp->fprintf("N=%-4d L=%-5d\n", num_nodes, num_links);
}


VOID
SCHEDULER::garbage_collect(LIST<HYP*>* hlist) {
  HYP *h;
  HYP *ph;
  if(Debug_flag.hyp)
    fprintf(stderr, "Garbage collecting: ");
  while(hlist->size > 0){
    h = hlist->pop();
    // delete any current HYPs that aren't extended, and chain backwards
    //  to delete dead-end paths
    while(h != NULL && h->get_n_next() == 0 && h != firsth) {
      ph = h->get_previous();
      if(Debug_flag.hyp)
	fprintf(stderr, "#");
      delete h;
      h = ph;
    }
  }
  if(Debug_flag.hyp)
    fprintf(stderr, "(%d left)\n", firsth->get_total_hyps());
}

VOID
SCHEDULER::hyp_chain_forward() {
  HYP* newfh;
  
  while(firsth->get_n_next() == 1) {
    firsth->online_display(dictionary->vocabulary, output_fp);
    if(state_decode_output) 
      firsth->online_state_display(dictionary->phonenames, firsth_start, 
				   frame_offset, frame_shift, state_dec_fp);
    if(ctm_output){
      if(Debug_flag.ndx) 
	ctm_fp->fprintf("**%.3f** ", 
		     (FLOAT)(start_time-firsth->get_time())*frame_shift);
      firsth->online_ctm_display(dictionary->vocabulary, firsth_start, 
				 frame_offset, frame_shift, ctm_fp, 
				 sentence_count, 'A');
    }
    if(srt_output)
      firsth->online_srt_display(dictionary->vocabulary, firsth_start, 
				 frame_offset, frame_shift, srt_fp);

    if(Verbose)
      firsth->online_log_display(dictionary->vocabulary, gstdout);	

    newfh = firsth->get_next()->pop();
    newfh->unset_previous();

    assert(outputs->size >= firsth->get_time() - firsth_start);
    for(INT tt = firsth_start; tt < firsth->get_time(); tt++){
      delete[] outputs->dequeue();
      lub->dequeue();
      def_lub_inc->dequeue();
    }

    firsth_start = firsth->get_time();
    delete firsth;
    firsth = newfh;
  }
  return;
}

	 
VOID
SCHEDULER::init_ndx(GENFILE* ndx_fp) {
  CHAR buf[BUFSIZ];
  CHAR bufcopy[BUFSIZ];
  CHAR* cptr;
  ndxq = new QUEUE<SECTION*>();
  // get the header
  if(ndx_fp->fgets(buf, BUFSIZ) == NULL)
    panic("Failed to read first line of ndx file!");
  ndx_hdr = strdup(buf);
  INT count = 1;
  FLOAT start, end;
  CHAR *typ, *tpc, *id;
  ndx_fp->fgets(buf, BUFSIZ);
  while(strncmp(buf, "</Episode>", 10)){
    strcpy(bufcopy, buf);
    if(ndx_fp->eof())
      panic("Premature eof of ndx file (no </Episode line)\n");

    cptr = strtok(buf, " \t"); // <Section
    if(strcmp(cptr, "<Section"))
      panic("Failed to read <Section in line %d of ndx file:\n %s [%s]", count, bufcopy, cptr);

    cptr = strtok(NULL, " \t"); // S_time=%d
    if(strncmp(cptr, "S_time=", 7))
      panic("Failed to read S_time in line %d of ndx file:\n %s [%s]", count, bufcopy, cptr);
    cptr += 7;
    start = atof(cptr);
    
    cptr = strtok(NULL, " \t"); // E_time=%d
    if(strncmp(cptr, "E_time=", 7))
      panic("Failed to read E_time in line %d of ndx file:\n %s [%s]", count, bufcopy, cptr);
    cptr += 7;
    end = atof(cptr);
    
    cptr = strtok(NULL, " \t"); // Type=%s
    if(strncmp(cptr, "Type=", 5))
      panic("Failed to read Type in line %d of ndx file:\n %s [%s]", count, bufcopy, cptr);
    cptr += 5;
    typ = strdup(cptr);

    cptr = strtok(NULL, " \t"); // Topic=%s
    if(strncmp(cptr, "Topic=", 6)){
      fprintf(stderr, "Failed to read Topic in line %d of ndx file:\n %s [%s]\n", count, bufcopy, cptr);
      tpc = strdup("");
    } else { 
      cptr += 6;
      tpc = strdup(cptr);
      cptr = strtok(NULL, " \t"); // ID=%s
    }

    if(strncmp(cptr, "ID=", 3))
      panic("Failed to read ID in line %d of ndx file:\n %s [%s]", count, bufcopy, cptr);
    cptr += 3;
    id = strdup(cptr);

    cptr = strtok(NULL, " \t"); // >
    if(strncmp(cptr, ">", 1))
      panic("Failed to read > in line %d of ndx file:\n %s [%s]", count, bufcopy, cptr);

    cptr = strtok(NULL, " \t");
    if(cptr != NULL)
      panic("Unexpected non-end of line %d in ndx file:\n %s [%s]", 
	    count, bufcopy, cptr);

    SECTION *s = new SECTION(start, end, frame_shift, typ, tpc, id);
    ndxq->enqueue(s);
    if(Debug_flag.ndx)
      s->info();
    delete typ; delete tpc; delete id;

    ndx_fp->fgets(buf, BUFSIZ);
  }
  ndx_ftr = strdup(buf); // </Episode>
  if(ndx_fp->fgets(buf, BUFSIZ) != NULL) {
    LONG nleft = (LONG) strlen(buf) - ndx_fp->ftell(); 
    ndx_fp->fseek(0L, SEEK_END);
    nleft += ndx_fp->ftell(); 
    fprintf(stderr, 
	    "%d bytes extra material after <\\Episode> in ndx file %s ignored\n", 
	    nleft, ndx_fp->get_name());
  }
}

