
#ifndef  include_acc
#define  include_acc

#include "inttypes.h"
#include "array.h"
#include "unit.h"

extern unsigned basemap[256]  ;

          /* accumulator: accumulates counts in bins
           * It uses EITHER acount OR hcount.
           * Thus all relevant code must be aware of that and check
           * which case is happening. In some cases
           * it is encapsulated, e.g. acc_bi_iterate.

           * In the calling code the acc struct itself is alloced on the
           * heap and initialized by acc_init.
          */
typedef struct          
{  rarray*acount        /* regular array count                             */
;  harray*hcount        /* harray count, who would have guessed            */
;  dim   n_mask         /* number of things not AGCTU                      */
;  dim   n_unsyl        /* number of missed sylamers                       */
;  dim   n_syl          /* number of passed sylamers                       */
;  dim   n_spaced       /* number of sylamers that fit, roughly            */
;  dim   n_unit         /* number of units in the window                   */
;  dim   n_r2           /* number of disallowed sylamers because of shift  */
;  dim   n_bases        /* total number of bases                           */
;  dim   k
;
}  acc   ;


void acc_null(acc* ac);
sylstatus acc_init(acc* ac, unsigned k, harray* ha);
void acc_reset(acc* ac);
void acc_free(acc* ac);
dim acc_get_sum(acc* ac);
void acc_update_diff(const acc* ac1, acc* ac2);
void mod_acc_init ( void );


sylstatus acc_addto
(  unsigned KK
,  unit* ut
,  acc* dst
,  int r2_check
)  ;

sylstatus acc_addto_gapped
(  unsigned KK
,  unit* ut
,  acc* dst
,  dim gap_min
,  dim gap_max
,  dim interval
,  acc* cache
)  ;

sylstatus acc_per_sequence
(  unsigned K
,  acc*  dst
,  acc*  src
,  unit* ut
,  dim   per_sequence_p
,  dim   cap
)  ;


unsigned sylid_from_buf(const char* buf, unsigned K, unsigned* m);
char* get_sylmer (int k, int sylid, char buf[16]);

void unit_filter
(  FILE* fp
,  const char* format
,  unit* ut
,  acc*  ac
,  dim count
)  ;

#endif

