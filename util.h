

#ifndef  include_util
#define  include_util

#include <stdlib.h>
#include <stdarg.h>

#include "inttypes.h"

#ifndef DEBUG
#define DEBUG  0
#endif

#define memdie die

void* mycalloc(size_t nmemb, size_t n_bytes);
void* mymalloc(size_t n_bytes);
void* myrealloc(void* mem, size_t n_bytes);

void tell(unsigned level,  const char* fmt,  ...);
void miaow(const char* fmt,  ...);        /* unsupressable warning/errors */
sylstatus arrr(const char* fmt,  ...)  ;  /* always returns SYL_FAIL   */
void* narrr(const char* fmt,  ...);       /* always returns NULL */
void die(const char* fmt,  ...)  ;

void myfname(char* buf, dim bufsize, const char* fmt, ...);

void frobotz(float g, float f, const char* s);

extern unsigned verbose_G        ;
extern unsigned debug_G          ;
extern unsigned n_alloc_exit_G   ;
extern double syl_db1_G          ;
extern unsigned char rmchar_G    ;

#define THEMIN(a,b)  (a < b ? a : b)
#define THEMAX(a,b)  (a > b ? a : b)

#define HASH_JUMP    149

   /*
    * other hash functions are generally abominable when word list is the full list
   */
#define hashwell(hid, mask)  (hid & mask)

#endif


