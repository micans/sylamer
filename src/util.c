
#include <stdio.h>

#include "util.h"
#include "inttypes.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#if DEBUG & DEBUG_MALLOC
#define REPORT(me, num) if(num>1024)miaow("%s %lu",me,(unsigned long)(num))
#else
#define REPORT(me, num)
#endif

unsigned verbose_G      =  1;
unsigned debug_G        =  0;
unsigned n_alloc_exit_G =  0;
double syl_db1_G        =  0.0;
unsigned char rmchar_G  = 'N';

static unsigned n_alloc_local = 0;

void tell
(  unsigned level
,  const char* fmt
,  ...
)
   {  va_list  args
   ;  if (verbose_G < level)
      return
   ;  va_start(args, fmt)
   ;  vfprintf(stderr, fmt, args)
   ;  fputc('\n', stderr)
   ;  va_end(args)
;  }


void* narrr
(  const char* fmt
,  ...
)
   {  va_list  args
   ;  va_start(args, fmt)
   ;  vfprintf(stderr, fmt, args)
   ;  fputc('\n', stderr)
   ;  va_end(args)
   ;  return NULL
;  }


sylstatus arrr
(  const char* fmt
,  ...
)
   {  va_list  args
   ;  va_start(args, fmt)
   ;  vfprintf(stderr, fmt, args)
   ;  fputc('\n', stderr)
   ;  va_end(args)
   ;  return SYL_FAIL
;  }


void miaow
(  const char* fmt
,  ...
)
   {  va_list  args
   ;  va_start(args, fmt)
   ;  vfprintf(stderr, fmt, args)
   ;  fputc('\n', stderr)
   ;  va_end(args)
;  }


void die
(  const char* fmt
,  ...
)
   {  va_list  args
   ;  va_start(args, fmt)
   ;  vfprintf(stderr, fmt, args)
   ;  fputc('\n', stderr)
   ;  va_end(args)
   ;  exit(1)
;  }


void* mycalloc(size_t nmemb, size_t n_bytes)
   {  void* mem
   ;  if (n_alloc_exit_G && n_alloc_local++ >= n_alloc_exit_G)
      mem = NULL
   ;  else
      mem = calloc(nmemb, n_bytes)

   ;  REPORT("mycalloc", nmemb * n_bytes)

   ;  if (nmemb && NULL == mem)
      miaow
      (  "___ memory shortage calloc asking for %lu bytes (calloc)"
      ,  (unsigned long) (n_bytes * nmemb)
      )
   ;  return mem
;  }


void* mymalloc(size_t n_bytes)
   {  void* mem
   ;  if (n_alloc_exit_G && n_alloc_local++ >= n_alloc_exit_G)
      mem = NULL
   ;  else
      mem = malloc(n_bytes)

   ;  REPORT("mymalloc", n_bytes)

   ;  if (n_bytes && NULL == mem)
      miaow
      (  "___ memory shortage malloc asking for %lu bytes (malloc)"
      ,  (unsigned long) n_bytes
      )
   ;  return mem
;  }


void* myrealloc(void* mem, size_t n_bytes)
   {  if (n_alloc_exit_G && n_alloc_local++ >= n_alloc_exit_G)
      mem = NULL
   ;  else
      mem = realloc(mem, n_bytes)
   ;  REPORT("myrealloc", n_bytes)
   ;  if (n_bytes > 0 && mem == NULL)
      miaow
      (  "___ memory shortage malloc asking for %lu bytes (realloc)"
      ,  (unsigned long) n_bytes
      )
   ;  return mem
;  }


void myfname
(  char* buf
,  dim bufsize
,  const char* fmt
,  ...
)
  {  va_list  args
  ;  va_start(args, fmt)

  ;  if (bufsize <= vsnprintf(buf, bufsize, fmt, args))
     miaow("cannot construct full file name -- it is truncated")
  ;  va_end(args)
; }


#if 0
void frobotz(float g, float f, const char* s) { fprintf(stderr, "%s frobotzes %g\n", s, f); }
#endif


