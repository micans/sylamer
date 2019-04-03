
#ifndef  include_interval
#define  include_interval

#include "unit.h"
#include "inttypes.h"

unsigned interval_prepare
(  collection* thecl
,  dim n_analyse
,  unsigned splitsize
,  unsigned binsize
,  int binbylength
)  ;


typedef struct
{  unsigned    i_hi_prev
;  unsigned    i_hi
;  unsigned    i_lo
;
}  spacer     ;


void interval_step
(  collection* thecl
,  dim n_analyse
,  spacer* sp
)  ;


#endif

