
#ifndef  include_sylio
#define  include_sylio

#include <stdio.h>
#include "inttypes.h"

#ifndef HAVE_ZLIB
#define HAVE_ZLIB 0
#endif

#if HAVE_ZLIB                 /* the INPUT abstraction carries some intricacies */
#include <zlib.h>
typedef gzFile INPUT;
#else
typedef FILE*  INPUT;
#endif

   /* if *status it failed;
    * if return value is NULL it reached end of file
   */
char* fread_line
(  INPUT fp
,  int* have_id
,  int flush_sequence
,  dim* linelen
,  dim* linecount
,  sylstatus* status
)  ;

   /* fixme document/check semantics
   */
char* fgets_line(INPUT fp, int* have_id, dim* ll, sylstatus* status);

void* myfopen (const char* fname, const char* mode);
void myfzclose(INPUT* fpp);
void myfclose(FILE** fpp);

#endif

