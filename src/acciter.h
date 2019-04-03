
#ifndef  include_acciter
#define  include_acciter

#include "inttypes.h"
#include "acc.h"
#include "array.h"

typedef struct
{  ofs      offset
;  dim      hid
;  dim      count1
;  dim      count2
;  float    val1
;  float    val2
;  float    exp
;
}  acc_iter ;

int acc_bi_iterate
(  acc_iter* ai
,  acc*  a1
,  acc*  a2
,  mkv*  m1
,  mkv*  m2
,  float* oracle
,  sylstatus* status
)  ;

int acc_iterate
(  acc_iter* ai
,  acc*  a1
,  sylstatus* status
)  ;

void acc_iter_reset(acc_iter* ai);

#endif

