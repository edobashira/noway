// -*- Mode: C++;  -*-
// File: decoder.cc
// Author: Steve Renals (s.renals@dcs.shef.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
// Copyright (C) Sheffield University, 1995
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  External wrapper main() function, 
//*             hooks into SCHEDULER::main_prog(argc, argv)
//*
//* CLASSES:
//* 
//* REQUIRES:  
//*
//* HISTORY:
//*  Jul 28 13:57 1994 (sjr): Added lattice headers and support for 
//*                           binary lattices
//*  Jul 27 12:20 1994 (sjr): Added timing information
//*  Jul 22 11:23 1994 (sjr): 0.4: implemented direct reading of lna files.
//*                           fixed up -tail_silence (changing name from
//*                           -trailing_silence).
//*  Jul 20 (sjr):  version 0.3;  added fast lookahead before 
//*                 activating new node (check for floor value for 
//*                 min duration frames ahead)
//*  Jul 20 (sjr):  release version 0.2;  includes support for lattice 
//*                 generation
//*  Jul 14 (sjr):  release version 0.1
//* Created: Fri Apr 22 14:09:20 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define DECLARE_GLOBAL_VARS
#include <memory>
#include <signal.h>
#include "all.h"
#include "NW_scheduler.h"

VOID
new_error() {
  panic("Operator new failed - out of memory\n");
}

#ifdef TIMING
VOID
interrupt_handler(INT) {
  if(Verbose){
    endt = times(&end_time_buffer);
    printf("\n\n");
    printf("Interrupt - Timing Information:\n");
    printf("  Elapsed real time = %lds\n", (endt-startt)/tick);
    printf("  User time = %lds\n", 
	   (end_time_buffer.tms_utime-start_time_buffer.tms_utime)/tick);
    printf("  System time = %lds\n", 
	   (end_time_buffer.tms_stime-start_time_buffer.tms_stime)/tick);
#if 0
#if defined(Linux) || defined(Solaris2)
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    printf("Memory Usage:\n");
    printf("  Max RSS = %ld, I RSS = %ld,  Data = %ld, Stack = %ld\n",
	   ru.ru_maxrss, ru.ru_ixrss, ru.ru_idrss, ru.ru_isrss);
    printf("  Swaps = %ld,  Page reclaims = %ld,  Page faults = %ld\n",
	   ru.ru_nswap, ru.ru_minflt, ru.ru_majflt);
#endif /* Linux Solaris2 */
#endif /* 0 */
    fflush(stdout);
  }  
  exit(1);
}
#endif

INT
main(INT argc, STRING argv[]) {
  SCHEDULER *sched = new SCHEDULER();
#ifdef TIMING
#if defined(_POSIX_VERSION)
  tick = sysconf(_SC_CLK_TCK);
#else
  tick = 100;
  fprintf(stderr, 
	 "Note assuming 100 ticks/s for timing purposes---this may be wrong!\n");
#endif /* _POSIX_VERSION */
  // if HP had implemented getrusage, I would change times to getrusage
  startt = times(&start_time_buffer);
  signal(SIGINT, interrupt_handler);
#endif
  std::set_new_handler(&new_error);
  gstdin = new GENFILE(stdin);
  gstdout = new GENFILE(stdout);
  gstderr = new GENFILE(stderr);
  sched->main_prog(argc, argv);
  exit(0);
}

