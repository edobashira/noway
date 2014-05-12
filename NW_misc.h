// -*- Mode: C++;  -*-
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "NW_type.h"

#ifndef MIN
#define MIN(X,Y) ((X)>(Y) ? (Y) : (X))
#endif
#ifndef MAX
#define MAX(X,Y) ((X)>(Y) ? (X) : (Y))
#endif
#define TWO_TO_THE(X) (1<<(int)(X))


// The following allows global variables to be declared in a header
// file without causing multiple definition errors in the linker
// and without having to manually maintain a separate file of
// definitions.
// The only file that defines DECLARE_GLOBAL_VARS is global.cc
// (and that is about all that it does).

#ifdef DECLARE_GLOBAL_VARS
#define GLOBAL_VAR(TYPE,NAME) TYPE NAME
#define GLOBAL_VAR_INIT(TYPE,NAME,INIT) TYPE NAME = INIT
#define GLOBAL_OBJ_INIT(CLASS,NAME,INIT) CLASS NAME INIT
#define EXTERN
#else
#define GLOBAL_VAR(TYPE,NAME) extern TYPE NAME
#define GLOBAL_VAR_INIT(TYPE,NAME,INIT) extern TYPE NAME
#define GLOBAL_OBJ_INIT(CLASS,NAME,INIT) extern CLASS NAME
#define EXTERN extern
#endif

// extern functions
class GENFILE;
extern VOID panic(const char*, ...);
extern INT logadd(INT loga, INT logb);
extern BOOL lnaread(FLOAT *probs, UCHAR *buffer, INT nout, GENFILE *fp);

const INT logadd_lookup_length = 0x8000;
EXTERN INT logadd_lookup[logadd_lookup_length];
const FLOAT log_scale = 8192.0;
const FLOAT one_over_log_scale = 0.00012207;
const DOUBLE l10 = 2.302585093;

inline DOUBLE
exp10(DOUBLE x) { return(exp(l10*x)); }

inline INT
logcode(FLOAT f) { return((INT)(log10(f) * log_scale + 0.5)); }

inline INT
logcodelog(FLOAT f) { return( (INT)(f * log_scale + 0.5)); }

inline FLOAT
logdecode(INT l) { return(exp10((FLOAT)l * one_over_log_scale)); }

#ifdef LNPROBS
 inline FLOAT
 logdecodelog(INT l) { return(((FLOAT)l * one_over_log_scale) * l10); }
#else
 inline FLOAT
 logdecodelog(INT l) { return((FLOAT)l * one_over_log_scale); }
 inline FLOAT
 logdecodeln(INT l) { return(((FLOAT)l * one_over_log_scale) * l10); }
#endif


const INT END_OF_SENTENCE = -1;

const INT PROB_NULL = 0x7f7f7f7f;
const USHORT LM_NULL = 0xffff;
const INT LM_OOV_MASK = 0x80000000;
const INT WORD_MASK   = ~0x80000000;

const UCHAR END_SENT_MASK = 0x80;

const INT ENTRY_STATE = 0;
const INT EXIT_STATE = 1;
const INT FIRST_STATE = 2;
const INT HEAD_IGNORE = 5;

GLOBAL_VAR(INT, Endian);
GLOBAL_VAR(BOOL, Verbose);
GLOBAL_VAR_INIT(BOOL, Aligning, FALSE);
GLOBAL_VAR_INIT(BOOL, Use_pron_priors, FALSE);

#ifndef DBL_MIN
#  define DBL_MIN   2.2250738585072014e-308   /* Min decimal value of a double */
#endif /* DBL_MIN */

#ifndef FLT_MIN 
#  define FLT_MIN   1.17549435e-38 /* Min decimal value of a float */
#endif /* FLT_MIN */

#ifndef DBL_DIG
#  define DBL_DIG   15     /* Digits of precision of a double */  
#endif /* DBL_DIG */

#ifndef DBL_MAX
#  define DBL_MAX   1.7976931348623157e+308   /* Max decimal value of a double */
#endif /* DBL_MAX */

#ifndef FLT_DIG
#  define FLT_DIG   6              /* Digits of precision of a float */  
#endif /* FLT_DIG */

#ifndef FLT_MAX
#  define FLT_MAX   3.40282347e+38 /* Max decimal value of a float */
#endif /* FLT_MAX */



class GENFILE {
private:
  FILE* fp;
  int pipe;
  char* name;
  size_t swapb(void* buf, size_t size, size_t nitems);  /* dpwe */
  int byteswap;       // flag to swap bytes of all accesses

public:
  GENFILE(const char* fname, const char* mode);

  GENFILE(FILE *f, int p = FALSE);

  ~GENFILE();

  FILE *get_fp() { return(fp); }

  char* get_name() { return name; }

  int fflush() { return(::fflush(fp)); }

  int fprintf(const char *format, ...);

  // Not all c libs have vfscanf() - so it's real hassle to define
  // fscanf here
  // the hack is to use fscanf(fil->get_fp(), ...) rather than fil->fscanf(...)

//  int fscanf(const char *format, ...) {
//    va_list args;
//    int ret;
//
//    va_start(args, format);
//    ret = vfscanf(fp, format, args);
//    va_end(args);
//    return(ret);
//  }

  int fgetc() {return(::fgetc(fp));}
  char* fgets(char *s, int n) {return(::fgets(s, n, fp));}
  int fputc(int c) {return(::fputc(c, fp));}
  int fputs(char *s) {return(::fputs(s, fp));}
  int ungetc(int c) {return(::ungetc(c, fp));}

  // dpwe added transparent byteswapping
  size_t fread(void *ptr, size_t size, size_t nobj){ 
    size_t n =::fread(ptr, size, nobj, fp);
    if(byteswap) swapb(ptr, size, n); 
    return n;
  }

  size_t fwrite(void *ptr, size_t size, size_t nobj){ 
    if(byteswap) swapb(ptr, size, nobj);
    size_t n=::fwrite(ptr, size, nobj, fp);
    if(byteswap) swapb(ptr, size, nobj); 
    return n;
  }

  //  size_t fread(void *ptr, size_t size, size_t nobj)
  //    { return(::fread(ptr, size, nobj, fp)); }
  //  size_t fwrite(void *ptr, size_t size, size_t nobj)
  //    { return(::fwrite(ptr, size, nobj, fp)); }

  int fseek(long int offset, int origin)
    { return(::fseek(fp, offset, origin)); }
  long int ftell() {return(::ftell(fp));}
  void rewind() {::rewind(fp);}
  //  int fgetpos(fpos_t *ptr) { return(::fgetpos(fp, ptr)); }
  //  int fsetpos(const fpos_t *ptr) { return(::fsetpos(fp, ptr)); }

  int eof() {return(feof(fp));}

  // this is needed when reading files written without swapping from
  // a littleendian machine
  void flip_byteswap()      
    { byteswap = !byteswap; }
    
};


GLOBAL_VAR(GENFILE*, gstdout);
GLOBAL_VAR(GENFILE*, gstdin);
GLOBAL_VAR(GENFILE*, gstderr);

