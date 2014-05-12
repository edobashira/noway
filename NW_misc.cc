// switch off assertions
// #define NDEBUG
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "NW_type.h"
#include "NW_misc.h"

VOID
panic(const char *format, ...)
{
    va_list args;

    va_start(args, format);
    fprintf(stderr, "Panic: ");
    (void) vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    (void) fflush(stderr);
    va_end(args);
    exit(1);
}

INT
logadd(INT a, INT b)
{
  static INT thresh = 0x7fff;
  INT diff = a-b > 0 ? a-b : b-a;
  INT m = MAX(a,b);
  if(diff > thresh){
    return(m);
  } else {
    return(m + logadd_lookup[diff]);
  }
}


/* Define statements */
const FLOAT LNPROB_FLOAT2INT = 24.0;

BOOL
lnaread(FLOAT *probs, UCHAR *buffer, INT nout, GENFILE *fp){
  //  FLOAT norm = 0.0;

  INT nread = fp->fread(buffer, sizeof(UCHAR), nout+1);

  if(fp->eof()) 
    return TRUE;

  if(nread != nout+1)
    panic("lna2float: Unexpected end of lna input\n");

  for(INT i = 0; i < nout; i++)
    probs[i] = exp( -(buffer[i+1] + 0.5) / LNPROB_FLOAT2INT);
  return((buffer[0] & END_SENT_MASK) != 0);
}

GENFILE::GENFILE(FILE *f, BOOL p) {
  fp = f;
  pipe = p;
  // check for byte order 
  INT i=1;
  CHAR* pc = (CHAR *)&i;
  byteswap = *pc;    // 0 for bigendian, 1 for littleendian 
  assert(sizeof(int) == 4);
  assert(sizeof(short) == 2);
}

GENFILE::GENFILE(STRING fname, STRING mode) {
  char buf[BUFSIZ];
  struct stat status;

  fp = NULL;
  pipe = FALSE;
  if(!strcmp(mode, "r")) {
    if(strlen(fname) > 4 && strcmp(&fname[strlen(fname) - 3], ".gz") == 0){
      if(stat(fname, &status) != 0) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
      sprintf(buf, "gzip --stdout --decompress --force %s", fname);
      if((fp = popen(buf, "r")) == NULL) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
      pipe = TRUE;
    } else if(strlen(fname) > 3 && strcmp(&fname[strlen(fname) - 2], ".Z") == 0) {
      if (stat(fname, &status) != 0) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
      sprintf(buf, "uncompress -c %s", fname);
      if ((fp = popen(buf, "r")) == NULL) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
      pipe = TRUE;
    } else {
      if ((fp = fopen(fname, "r")) == NULL) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
    }
  } else if (!strcmp(mode, "w")) {
    if (strlen(fname) > 4 && strcmp(&fname[strlen(fname) - 3], ".gz") == 0 ){
      sprintf(buf, "gzip --stdout --force  > %s", fname);
      if ((fp = popen(buf, "w")) == NULL) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
      pipe = TRUE;
    } else if(strlen(fname) > 3 && strcmp(&fname[strlen(fname) - 2], ".Z") == 0) {
      sprintf(buf, "compress -c > %s", fname);
      if ((fp = popen(buf, "w")) == NULL) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
      pipe = TRUE;
    } else {
      if ((fp = fopen(fname, "w")) == NULL) {
	perror("GENFILE::GENFILE()");
	panic("cannot open '%s'\n", fname);
      }
    }
  } else {
    if ((strlen(fname) > 3 && strcmp(&fname[strlen(fname) - 3], ".gz") == 0) ||
	(strlen(fname) > 2 && strcmp(&fname[strlen(fname) - 2], ".Z") == 0)) 
      ::fprintf(stderr, "WARNING: Cannot open named pipe using mode %s (only r or w)\n", mode);
    if ((fp = fopen(fname, mode)) == NULL) {
      perror("GENFILE::GENFILE(STRING, STRING)");
      panic("cannot open '%s'\n", fname);
    }
  }

  name = strdup(fname);

  // check for byte order 
  INT i=1;
  CHAR* pc = (CHAR *)&i;
  byteswap = *pc;    // 0 for bigendian, 1 for littleendian 
  assert(sizeof(int) == 4);
  assert(sizeof(short) == 2);
}


GENFILE::~GENFILE() {
  if(!(fp == stdin || fp == stderr || fp == stdout || fp == NULL)){
    if(pipe)
      pclose(fp);
    else
      fclose(fp);
  }
  delete name;
}


INT 
GENFILE::fprintf(const char *format, ...) {
  va_list args;
  INT ret;
  
  va_start(args, format);
  ret = vfprintf(fp, format, args);
  va_end(args);
  return(ret);
}

size_t GENFILE::swapb(void *buf, size_t size, size_t nitems) {
  // swap pairs or quads of bytes in a buffer 
  if(size==2) {
    SHORT *sbuf = (SHORT *)buf;
    INT x;
    // assert(sizeof(short)==2);
    while(nitems--) {
      x = *sbuf;
      *sbuf++ = ((x & 0xFF) << 8) + ((x >> 8) & 0xFF);
    }
  } else if(size==4) {
    INT *ibuf = (INT *)buf;
    INT x;
    // assert(sizeof(int)==4);
    while(nitems--) {
      x = *ibuf;
      *ibuf++ = ((x & 0xFFL) << 24L) + ((x & 0xFF00L) << 8L) \
	+ ((x & 0xFF0000L) >> 8L) + ((x >> 24L) & 0xFF);
    }
  }
  // ignore other sizes 
  return nitems;
}
