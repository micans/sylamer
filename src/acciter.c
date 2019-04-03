

#include "inttypes.h"
#include "util.h"
#include "acciter.h"
#include "acc.h"
#include "array.h"


      /* this iterator abstracts over how the counts are stored,
       * either in an array or in a harray.
       * m1 and m2 are either both NULL or both valid.
       *
       * Oracle is sorted in the same order as the word list,
       * (if the oracle is present). This is a somewhat loose
       * and unsatisfactory dependency. (but note that they
       * are read in simultaneously).
      */
int acc_bi_iterate
(  acc_iter* ai
,  acc*  a1
,  acc*  a2
,  mkv*  m1
,  mkv*  m2
,  float* oracle
,  sylstatus* mystatus
)
   {  unsigned hid
   ;  unsigned hhh
   ;  dim count
   ;  dim os  

   ;  ai->offset
      =  array_step
         (  ai->offset
         ,  a1->hcount
         ,  a1->acount
         ,  &hid
         ,  &hhh
         ,  &count
         ,  NULL
         ,  mystatus
         )
   ;  if (ai->offset < 0)
      return 0

   ;  os = ai->offset-1
   ;  ai->hid = hid

   ;  if (a1->acount)
      {  ai->count1  =  count
      ;  ai->count2  =  a2->acount->cells[os].num
      ;  if (m1)
            ai->val1 =  m1->rval->cells[os].val
         ,  ai->val2 =  m2->rval->cells[os].val
      ;  else
            ai->val1 =  0.0
         ,  ai->val2 =  0.0
   ;  }
      else if (a1->hcount)
      {  hcell* cells1 = a1->hcount->cells
      ;  hcell* cells2 = a2->hcount->cells
      ;  hcell* cells3 = m1 ? m1->hval->cells : NULL
      ;  hcell* cells4 = m1 ? m2->hval->cells : NULL

      ;  if
         (  cells1[hhh].hid != hid
         || cells2[hhh].hid != hid
         || (m1 && cells3[hhh].hid != hid)
         || (m1 && cells4[hhh].hid != hid)
         || cells1[hhh].pl.num != count
         )
         {  *mystatus = SYL_FAIL  
         ;  miaow
            (  "hash bg/window conflict [hash %lu, hids %lu,%lu, markov %s]"
            ,  hhh, cells1[hhh].hid, cells2[hhh].hid
            ,  m1 ? "not shown" : "not used"
            )
         ;  return 0
      ;  }

         ai->count1  =  cells1[hhh].pl.num
      ;  ai->count2  =  cells2[hhh].pl.num
      ;  if (m1)
            ai->val1 =  cells3[hhh].pl.val
         ,  ai->val2 =  cells4[hhh].pl.val
      ;  else
            ai->val1 =  0.0
         ,  ai->val2 =  0.0
      ;  if (oracle)
         ai->exp = oracle[os]
      ;  else
         ai->exp = 0.0
   ;  }
      return 1
;  }


      /* iterate over a single acc instance */
int acc_iterate
(  acc_iter* ai
,  acc*  a1
,  sylstatus *mystatus
)
   {  unsigned hid
   ;  unsigned hhh
   ;  dim count
   ;  dim os  

   ;  ai->offset
      =  array_step
         (  ai->offset
         ,  a1->hcount
         ,  a1->acount
         ,  &hid
         ,  &hhh
         ,  &count
         ,  NULL
         ,  mystatus
         )
   ;  if (ai->offset < 0)
      return 0

   ;  os = ai->offset-1
   ;  ai->hid = hid

   ;  if (a1->acount)
      ai->count1  =  count
   ;  else if (a1->hcount)
      ai->count1  =  a1->hcount->cells[hhh].pl.num

   ;  return 1
;  }


void acc_iter_reset(acc_iter* ai)
   {  ai->offset  =  0
   ;  ai->hid     =  0
   ;  ai->count1  =  0
   ;  ai->count2  =  0
   ;  ai->val1    =  0.0
   ;  ai->val2    =  0.0
   ;  ai->exp     =  0.0
;  }


