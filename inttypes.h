

#ifndef  include_inttypes
#define  include_inttypes

#include <limits.h>
#include <sys/types.h>

/* Some vital derived constants
 * LSIZE    (language size) is the number of words of length k over 4-ary alphabet.
 * LOMEGA   is the last (all bits one) word. It is used as a mask in places.
 * MERMASK  is used to check validity of the current sylamer. It uses a single bit per position.
*/

#define LSIZE(k)     (1 << (2*(k)))    /* language size, 4096, 16384          */
#define LOMEGA(k)    (LSIZE(k)-1)      /* last (all-one) word, 4095, 16383,   */
#define MERMASK(k)   ((1 << (k))-1)    /* 6->63 and 7->127 (K-bits-all-one)   */


enum
{  VERBOSE_ONCE = 1
,  VERBOSE_BIN  = 2
,  VERBOSE_MAX  = 3   
}  ;

typedef unsigned long dim;
typedef long ofs;

#ifndef ulong
#  define ulong unsigned long
#endif

#define DIM_MAX ULONG_MAX

typedef unsigned sylstatus;      /* zero indicates success */

enum
{  SYL_SUCCESS =  0u
,  SYL_FAIL    =  1u
,  SYL_DONE
,  SYL_IGNORE
}  ;

#define REP_OVER  1
#define REP_UNDER 2

#define STOP_BG   1
#define STOP_READ 2

#define DEBUG_READ 1
#define DEBUG_MALLOC 2

#endif


