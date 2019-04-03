

#include "util.h"
#include "acc.h"
#include "array.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void rarray_free(rarray* ra)
   {  free(ra->cells)
   ;  free(ra)
;  }


void harray_free(harray* ha)
   {  free(ha->cells)
   ;  free(ha)
;  }


harray* harray_clone (harray* ha)
   {  dim i
   ;  harray* hb = mymalloc(sizeof hb[0])
   ;  if (hb)
      hb->cells = mymalloc((ha->hmask+1) * sizeof hb->cells[0])

   ;  if (!hb || !hb->cells)      /* we leak hb, but it's order sizeof(void*) */
      return NULL

   ;  hb->hmask = ha->hmask
   ;  hb->words = ha->words
   ;  hb->n_words= ha->n_words

   ;  for (i=0;i<=ha->hmask;i++)
      {  hb->cells[i].pl.num = 0
      ;  hb->cells[i].hid = ha->cells[i].hid
   ;  }
      return hb
;  }



void harray_reset (harray* ha)
   {  dim i
   ;  for (i=0;i<=ha->hmask;i++)
      ha->cells[i].pl.num = 0
;  }


void rarray_reset (rarray* ra)
   {  memset(ra->cells, 0, ra->n_cells * sizeof ra->cells[0])
;  }


rarray* rarray_new(unsigned k)
   {  rarray* ra  =  mymalloc(sizeof ra[0])
   ;  if (ra)
      ra->cells = mycalloc(LSIZE(k), sizeof ra->cells[0])

   ;  if (!ra || !ra->cells)        /* we might near-infinitesimally leak ra */
      return NULL
   ;  ra->n_cells =  LSIZE(k)
   ;  return ra
;  }


void mkv_null
(  mkv* markov
)
   {  markov->hval =  NULL
   ;  markov->rval =  NULL
;  }


sylstatus mkv_init
(  mkv* markov
,  unsigned k
,  harray* ha
)
   {  if (ha)
         markov->hval =  harray_clone(ha)
      ,  markov->rval =  NULL
   ;  else
         markov->hval =  NULL
      ,  markov->rval =  rarray_new(k)
   ;  return markov->hval || markov->rval ? SYL_SUCCESS : SYL_FAIL
;  }


void mkv_free
(  mkv* markov
)
   {  if (markov->rval)
      rarray_free(markov->rval)
   ;  if (markov->hval)
      harray_free(markov->hval)
;  }


harray* harray_new (wlist* wl, unsigned mask, unsigned K)
   {  dim i
   ;  harray* ha  =  mymalloc(sizeof ha[0])
   ;  hcell* cp   =  NULL
   ;  unsigned h_miss = 0
   ;  char buf[16]

   ;  if (ha)
      cp =  ha->cells
         =  mymalloc((mask + 1) * sizeof ha->cells[0])

   ;  if (!ha || !ha->cells)
      return NULL

   ;  ha->hmask = mask

   ;  for (i=0; i<=mask;i++)
         cp[i].hid = -1u         /* NOTE this is TTTT,TTTT,TTTT,TTTT (assuming 32 bits) */
      ,  cp[i].pl.num = 0

                                 /* fixme: dependency with code in acc.c */
   ;  for (i=0; i<wl->n_words;i++)
      {  unsigned hhh = hashwell(wl->words[i], mask)
      ;  dim jj = 0
      ;  while(cp[hhh].hid != -1u && cp[hhh].hid != wl->words[i])
         {  hhh += HASH_JUMP    /* not a factor of 2^k   */
         ;  hhh &= mask
         ;  if (jj++ > mask)
            {  miaow("hash error in harray_new")
            ;  free(cp)
            ;  free(ha)
            ;  return NULL
         ;  }
            h_miss++
      ;  }
         cp[hhh].pl.num = 0
      ;  if (cp[hhh].hid == wl->words[i])
         miaow("___[warning] word %s present more than once", get_sylmer(K, wl->words[i], buf))
      ;  cp[hhh].hid = wl->words[i]
   ;  }

      ha->words = wl->words
   ;  ha->n_words = wl->n_words
   ;  tell
      (  VERBOSE_ONCE
      ,  "[hash] %u words, word mask %u, hash misses %u"
      ,  (unsigned) wl->n_words
      ,  mask
      ,  (unsigned) h_miss
      )
   ;  return ha
;  }


ofs array_step
(  ofs         offset
,  harray*     ha
,  rarray*     ra
,  unsigned*   hidp
,  unsigned*   hhhp        /* the pivotal hash index, if applicable */
,  dim*        countp
,  float**     valpp
,  sylstatus   *mystatus
)
   {  const unsigned* words = NULL
   ;  dim n_words = 0

   ;  if (ha)
         words =  ha->words
      ,  n_words = ha->n_words
      
   ;  if (valpp)
      *valpp = NULL

   ;  if (hhhp)
      *hhhp = 0

   ;  if ((words && offset >= n_words) || (ra && offset >= ra->n_cells))
      return -1

   ;  if (words)
      {  hcell* cells = ha->cells
      ;  unsigned hid, hhh, hmask
      ;  dim jj = 0

      ;  hid = words[offset]
      ;  hmask = ha->hmask
      ;  hhh = hashwell(hid, hmask)

      ;  while(cells[hhh].hid != -1u && cells[hhh].hid != hid)
         {  hhh += HASH_JUMP
         ;  hhh &= hmask
         ;  if (jj++ > hmask)
            {  *mystatus = SYL_FAIL
            ;  miaow("hash error in array_step")
            ;  return -1
         ;  }
      ;  }
         if (cells[hhh].hid != hid)
         {  *mystatus = SYL_FAIL
         ;  miaow
            (  "hash markov conflict [hash %lu, hids %lu/%lu]"
            ,  hhh
            ,  cells[hhh].hid
            ,  hid
            )
         ;  return -1
      ;  }

         *hidp = hid
      ;  if (hhhp)
         *hhhp = hhh

      ;  if (valpp)
         *valpp = &(cells[hhh].pl.val)
      ;  else if (countp)                 /* exclusive, it's a union */
         *countp = cells[hhh].pl.num
   ;  }
      else
      {  *hidp = offset          /* beware, signed in unsigned */
      ;  if (valpp)
         *valpp = &(ra->cells[offset].val)
      ;  else if (countp)
         *countp = ra->cells[offset].num
   ;  }
      return offset + 1
;  }



ofs mkv_iterate
(  ofs         offset      /* for either markov->hval or markov->rval */
,  mkv*        markov
,  unsigned*   hidp        /* store the current word                  */
,  float**     valpp       /* make available for writing              */
,  sylstatus   *mystatus
)
   {  ofs o
      =  array_step
         (  offset
         ,  markov->hval
         ,  markov->rval
         ,  hidp
         ,  NULL     /* not interested in hash index */
         ,  NULL     /* there is no count in this case */
         ,  valpp    /* but rather storage for our epxected frequency */
         ,  mystatus
         )
   ;  return o
;  }


               /* dead code. dangling on the bick bucket rim */
#if 0
int hcell_cmp(const void* p1, const void* p2)
   {  const hcell*h1 = p1, *h2 = p2
   ;  return h1->hid < h2->hid ? -1 : h1->hid == h2->hid ? 0 : 1
;  }

void harray_sort(harray* ha)
   {  qsort(ha->cells, ha->hmask+1, sizeof ha->cells[0], hcell_cmp)
;  }
#endif


