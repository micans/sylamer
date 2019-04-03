

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "acc.h"
#include "inttypes.h"
#include "array.h"
#include "unit.h"
#include "util.h"


/* TODO
   -  factor out commonality in hash code
   -  ship harray rather than cells + mask
*/

unsigned basemap[256] = { 0 };

void mod_acc_init ( void )
   {  dim i
   ;  for (i=0;i<256;i++)
      basemap[i] = 5
   ;  basemap['A'] = 0
   ;  basemap['C'] = 1
   ;  basemap['G'] = 2
   ;  basemap['T'] = 3
   ;  basemap['U'] = 3     /* everything < 4 is coding */
   ;  basemap['N'] = 4     /* masked */
   ;  basemap['X'] = 4     /* masked */

   ;  basemap['\n'] = 8    /* bit 8 set */
   ;  basemap['\r'] = 9    /* bit 8 set */
;  }


void acc_null
(  acc* ac
)
   {  ac->acount = NULL
   ;  ac->hcount = NULL
   ;  ac->n_unsyl = 0
   ;  ac->n_mask = 0
   ;  ac->n_syl = 0
   ;  ac->n_spaced = 0
   ;  ac->n_unit = 0
   ;  ac->n_r2 = 0
   ;  ac->n_bases = 0
   ;  ac->k = 0
;  }


void acc_update_diff
(  const acc* ac1
,  acc* ac2
)
   {  ac2->n_unsyl   =  ac1->n_unsyl  - ac2->n_unsyl
   ;  ac2->n_mask    =  ac1->n_mask   - ac2->n_mask
   ;  ac2->n_syl     =  ac1->n_syl    - ac2->n_syl
   ;  ac2->n_spaced  =  ac1->n_spaced - ac2->n_spaced
   ;  ac2->n_unit    =  ac1->n_unit   - ac2->n_unit
   ;  ac2->n_r2      =  ac1->n_r2     - ac2->n_r2
   ;  ac2->n_bases   =  ac1->n_bases  - ac2->n_bases
;  }


void acc_free
(  acc* ac
)
   {  if (ac->acount)
      rarray_free(ac->acount)
   ;  if (ac->hcount)
      harray_free(ac->hcount)
;  }


               /* when stepping bins rather than growing, we reset */
void acc_reset
(  acc* ac
)
   {  if (ac->acount)
      rarray_reset(ac->acount)
   ;  if (ac->hcount)
      harray_reset(ac->hcount)
   ;  ac->n_unsyl = 0
   ;  ac->n_mask = 0
   ;  ac->n_syl = 0
   ;  ac->n_spaced = 0
   ;  ac->n_unit = 0
   ;  ac->n_r2 = 0
   ;  ac->n_bases = 0
;  }


sylstatus acc_init
(  acc* ac
,  unsigned k
,  harray* ha
)
   {  if (ha)
         ac->hcount = harray_clone(ha)
      ,  ac->acount = NULL
   ;  else
         ac->acount = rarray_new(k)
      ,  ac->hcount = NULL

   ;  ac->n_unsyl =  0
   ;  ac->n_mask  =  0
   ;  ac->n_syl   =  0
   ;  ac->n_spaced=  0
   ;  ac->n_unit  =  0
   ;  ac->n_r2    =  0
   ;  ac->n_bases =  0
   ;  ac->k       =  k

   ;  return ac->hcount || ac->acount ? SYL_SUCCESS : SYL_FAIL
;  }


dim acc_get_sum(acc* ac)
   {  dim i, sum = 0
   ;  if (ac->acount)
      for (i=0;i<ac->acount->n_cells;i++)
      sum += ac->acount->cells[i].num
   ;  else
      for (i=0;i<=ac->hcount->hmask;i++)
      sum += ac->hcount->cells[i].pl.num
   ;  return sum
;  }


sylstatus acc_addto
(  unsigned K
,  unit* ut
,  acc* dst
,  int r2_check
)
   {  dim i, n_mask = 0, n_unsyl = 0, n_syl = 0, n_spaced = 0, n_r2 = 0
   ;  const unsigned char* bases =  ut->bases
   ;  const dim n_bases =  ut->n_bases
   ;  rcell* dst_arry   =  dst->acount ? dst->acount->cells : NULL
   ;  harray* dst_hash  =  dst->hcount

   ;  unsigned hid      =  0
   ;  unsigned mask     =  0
   ;  unsigned stretchlength = 0

   ;  unsigned cache[2] = { -1, -1 }

   ;  if (!dst_arry && !dst_hash)
      return arrr("programmer bd in acc_addto")

   ;  if (!K)
      return SYL_SUCCESS         /* recentlychanged when relaxing and allowing m = 1 */

   ;  if (debug_G && bases)
      {  int slen = strlen((char*) bases)
      ;  if (slen != n_bases)
         return arrr
         (  "sanity check on input failed for sequence <%s>: "
            "measured/passed %d vs %d"
         ,  ut->sid
         ,  slen
         ,  n_bases
         )
   ;  }

      if (!n_bases)
      tell(VERBOSE_BIN, "warning: empty sequence for <%s>", ut->sid)

   ;  if (n_bases < K)
      return SYL_SUCCESS

                  /* initialize the mapping: compute first sylamer map.
                   * FIXME: simply use K-1, forget about adding the count,
                   * similar to dyad code
                  */
   ;  for (i=0;i<K;i++)
      {  unsigned b = basemap[bases[i]]
      ;  hid |= b << 2*(K-1-i)
      ;  mask = mask << 1
      ;  if (b >= 4)
         {  mask |= 1
         ;  n_mask++
         ;  stretchlength = 0
      ;  }
         else
         stretchlength++
   ;  }


      if (mask & MERMASK(K))
      n_unsyl++
   ;  else
      {  if (dst_arry)
         dst_arry[hid].num++
      ;  else
         HINC_IF_PRESENT(dst_hash, hid)
      ;  cache[0] = hid
      ;  n_syl++
   ;  }

      for (i=K;i<n_bases;i++)
      {  unsigned b = basemap[bases[i]]

      ;  mask = mask << 1
      ;  hid = ((hid << 2) | b) & LOMEGA(K)

      ;  if (b >= 4)                /* b is masked/invalid */
         {  mask |= 1
         ;  n_mask++
         ;  n_spaced += stretchlength/K
         ;  stretchlength = 0
      ;  }
         else
         stretchlength++

      ;  if (mask & MERMASK(K))
         n_unsyl++
                        /* mq generalize this to arbitrary lengths */
      ;  else
         {  if (r2_check && (hid == cache[0] || hid == cache[1]))
            {  n_r2++            /* mqr2 */
            ;  continue
         ;  }
            cache[i&1] = hid
         ;  n_syl++
         ;  if (dst_arry)
            dst_arry[hid].num++
         ;  else
            HINC_IF_PRESENT(dst_hash, hid)
      ;  }
      }
      dst->n_mask  += n_mask
   ;  dst->n_unsyl += n_unsyl
   ;  dst->n_syl   += n_syl
   ;  dst->n_spaced   += n_spaced + stretchlength / K
   ;  dst->n_unit++
   ;  dst->n_r2    += n_r2
   ;  dst->n_bases += n_bases
   ;  return SYL_SUCCESS
;  }


   /* This one adds per-unit count to dst, looking at those counts in src
    * that exceed per_sequence_p.
    * It does so by looping over the bases in the unit.
    * While looping, it clears the counts in src. This is for two reasons
    * important: First, src has to be reset. Second, the count has
    * to be only checked once. The second time we inspect count (because
    * the corresponding k-mer was seen more than once in ut) we will
    * do a no-op because the count was reset -- and that's what we want.
    *
    * Some code duplication with acc_addto, but unifying them
    * would result in horrendous spagetti. Perhaps macrolizing
    * or functionating some of the shared blocks can be considered.
   */
sylstatus acc_per_sequence
(  unsigned K
,  acc*  dst
,  acc*  src
,  unit* ut
,  dim   per_sequence_p
,  dim   cap
)
   {  dim i, n_mask = 0
   ;  const unsigned char* bases =  ut->bases
   ;  const dim n_bases =  ut->n_bases

   ;  rcell* src_arry   =  src->acount ? src->acount->cells : NULL
   ;  harray* src_hash  =  src->hcount

   ;  rcell* dst_arry   =  dst->acount ? dst->acount->cells : NULL
   ;  harray* dst_hash  =  dst->hcount

   ;  unsigned hid      =  0
   ;  unsigned mask     =  0
   ;  unsigned count    =  0

   ;  dst->n_unit  += 1

   ;  if (n_bases < K)
      return SYL_SUCCESS

                  /* initialize the mapping: compute first sylamer map */
   ;  for (i=0;i<K;i++)
      {  unsigned b = basemap[bases[i]]
      ;  hid |= b << 2*(K-1-i)
      ;  mask = mask << 1
      ;  if (b >= 4)
         {  mask |= 1
         ;  n_mask++
      ;  }
      }

      if (!(mask & MERMASK(K))) 
      {  if (src_arry)
         {  count = src_arry[hid].num
         ;  src_arry[hid].num = 0
      ;  }
         else
         HZERO_IF_PRESENT(count, src_hash, hid)

      ;  {  dim inc = 0
         ;  if (per_sequence_p && count >= per_sequence_p)
            inc = 1
         ;  else if (cap)
            {  inc = count
            ;  if (count > cap)
               count = cap
         ;  }

            if (dst_arry)
            dst_arry[hid].num += inc
         ;  else
            HADD_IF_PRESENT(dst_hash, hid, inc)
      ;  }
      }

      for (i=K;i<n_bases;i++)
      {  unsigned b = basemap[bases[i]]

      ;  mask = mask << 1
      ;  hid = ((hid << 2) | b) & LOMEGA(K)

      ;  if (b >= 4)                /* b is masked/invalid */
         {  mask |= 1
         ;  n_mask++
      ;  }

         if (!(mask & MERMASK(K)))
         {  if (src_arry)
            {  count = src_arry[hid].num
            ;  src_arry[hid].num = 0
         ;  }
            else
            HZERO_IF_PRESENT(count, src_hash, hid)

         ;  {  dim inc = 0
            ;  if (per_sequence_p && count >= per_sequence_p)
               inc = 1
            ;  else if (cap)
               {  inc = count
               ;  if (inc > cap)
                  inc = cap
            ;  }

               if (dst_arry)
               dst_arry[hid].num += inc
            ;  else
               HADD_IF_PRESENT(dst_hash, hid, inc)
         ;  }
         }
if (0 && hid == 15 && count > 0)fprintf(stderr, "15 -> %d %d\n", (int) count, (int) dst_arry[hid].num);
      }
      dst->n_mask  += src->n_mask
   ;  dst->n_unsyl += src->n_unsyl
   ;  dst->n_spaced+= src->n_spaced
   ;  dst->n_r2    += src->n_r2
   ;  dst->n_bases += src->n_bases
   ;  dst->n_syl   += src->n_syl

   ;  src->n_mask   = 0
   ;  src->n_unsyl  = 0
   ;  src->n_syl    = 0
   ;  src->n_spaced = 0
   ;  src->n_r2     = 0
   ;  src->n_bases  = 0
   ;  src->n_unit   = 0
   ;  return SYL_SUCCESS
;  }


   /* Get hash from a word
   */
unsigned sylid_from_buf(const char* buf, unsigned K, unsigned* m)
   {  unsigned hid = 0
   ;  unsigned mask = 0
   ;  dim i
   ;  for (i=0;i<K;i++)
      {  unsigned b = basemap[(unsigned char) buf[i]]
      ;  hid |= b << 2*(K-1-i)
      ;  mask = mask << 1
      ;  if (b >= 4)
         mask |= 1
   ;  }
      *m = mask
   ;  return hid
;  }

static char getbase[4] = { 'A', 'C', 'G', 'T' };

         /* Get the word that corresponds to a hash
         */
char* get_sylmer
(  int k
,  int sylid
,  char buf[16]
)
   {  int i

   ;  if (sylid < 0 || sylid > LOMEGA(k))
      {  if(0)miaow("get_sylmer[%d] has corrupted input <%d>", k, sylid)
      ;  for (i=0;i<k;i++)
         buf[i] = 'X'
   ;  }
      else
      for (i=0;i<k;i++)
      buf[i] = getbase[(sylid >> (2*(k-i-1))) & 3]

   ;  buf[k] = '\0'
   ;  return buf
;  }


static void put_revcompl
(  const unsigned char* bases
,  dim n_bases
,  FILE* fp
,  int format
)
   {  dim j
   ;  for (j=0;j<n_bases;j++)
      {  unsigned char c = bases[n_bases-j-1]
      ;  unsigned char t =      c == 'T'             /* fixme; lookup table? */
                     ?  'A'
                     :     c == 'A'
                        ?  'T'
                        :     c == 'G'
                           ?  'C'
                           :     c == 'C'
                              ?  'G'
                              :  c
      ;  if (format && j && !(j%60)) fputc('\n', fp)
      ;  fputc(t, fp)
   ;  }
   }


void unit_filter_fragment
(  FILE* fp
,  unsigned char* bases
,  dim n_bases
)
   {  dim n_runs     =  0
   ;  dim max_size   =  0
   ;  dim min_size   =  DIM_MAX
   ;  dim n_base     =  0
   ;  dim nn         = 0
   ;  int run_length =  0
   ;  dim j

   ;  for (j=0; j<n_bases; j++)
      {  if (basemap[bases[j]] > 3)
         {  if (run_length)
            {  if (run_length > max_size)  max_size = run_length 
            ;  if (run_length < min_size)  min_size = run_length 
            ;  n_runs++
            ;  run_length = 0
         ;  }
            nn++
      ;  }
         else
            run_length++
         ,  n_base++
   ;  }
      fprintf
      (  fp
      ,  "nfrag=%lu,max=%lu,min=%lu,nbase=%lu,nn=%lu"
      ,  (ulong) n_runs
      ,  (ulong) max_size
      ,  (ulong) min_size
      ,  (ulong) n_base
      ,  (ulong) nn
      )
;  }


void unit_filter
(  FILE* fp
,  const char* format
,  unit* ut
,  acc*  ac
,  dim count
)
   {  const char* f = format
   ;  const char* z = format + strlen(format)
   ;  int escape = 0
   ;  while (f < z)
      {  unsigned int c = (unsigned char) f[0]
      ;  if (!escape)
         {  if (c == '%')
            escape = 1
         ;  else
            fputc(c, fp)
      ;  }
         else
         {  switch(c)
            {  case 'S': fputs((const char*) ut->bases, fp); break
            ;  case 'T': put_revcompl(ut->bases, ut->n_bases, fp, 0); break
            ;  case 'G': put_revcompl(ut->bases, ut->n_bases, fp, 1); break
            ;  case 'A': fputs(ut->annot, fp); break
            ;  case 'I': fputs(ut->sid, fp); break
            ;  case 'C': fprintf(fp, "%lu", (unsigned long) count); break
            ;  case 'L': fprintf(fp, "%lu", (unsigned long) ut->n_bases); break
            ;  case 'X': unit_filter_fragment(fp, ut->bases, ut->n_bases); break
            ;  case 'F': { dim j
                         ; for(j=0;j<ut->n_bases;j++)
                           {  if (j && !(j%60)) fputc('\n', fp)
                           ;  fputc(ut->bases[j], fp)
                         ; }  
                         } break
            ;  case 'R': fprintf(fp, "%lu", (unsigned long) ut->rank); break
            ;  case 'n': fputc('\n', fp); break
            ;  case 't': fputc('\t', fp); break
            ;  case '%': fputc('%', fp); break
            ;  default : fprintf(fp, "%c", c); break
         ;  }
            escape = 0
      ;  }
         f++
   ;  }
   }


/* fixme: repeat checking as implemented does not make sufficient sense.
 * e.g.
AGCGGGTGCACACACGTCAC
      TGC   CAC
we skip this because CAC is a 2-repeat, yet we've not seen TGC-CAC
before because the distance is too small.
Similarly, if we skip a 2-repeat left dyad, we may set the count
for certain combinations to zero.
So for now, repeat skipping has been disabled. Perhaps it needs
to be implemented at the conjoined-dyad level, although a little
consideration shows it to be a non-trivial problem.

It's hard to just do with a linear list of previous stored words.
We'll have to store the last average offset of the two dyads,
in a separate array.

Update 10-208. Implemented simple repeat check on both hid1 and hid1.
Update 10-209. Removed said simple repeat check.

*/


sylstatus acc_addto_gapped
(  unsigned KK
,  unit* ut
,  acc* dst
,  dim gap_min
,  dim gap_max
,  dim interval
,  acc* cache
)
   {  dim l, r
   ;  const unsigned char* bases =  ut->bases
   ;  const dim n_bases =  ut->n_bases
   ;  rcell* dst_arry   =  dst->acount ? dst->acount->cells : NULL

   ;  unsigned n_unsyl  =  0, n_syl = 0, n_r2 = 0, n_mask = 0
   ;  const char* space =  "                                                                          "

   ;  unsigned hid1 = 0, mask1 = 0
   ;  unsigned hid2 = 0, mask2 = 0

   ;  unsigned b1, b2
   ;  unsigned K = KK/2

   ;  char buf1[16], buf2[16]

   ;  if (!KK || (KK&1))
      return arrr("need positive even K")

   ;  if (!dst_arry && !dst->hcount)
      return arrr("programmer bd in acc_addto_gapped")

   ;  if (!n_bases)
      tell(VERBOSE_BIN, "warning: empty sequence for <%s>", ut->sid)

   ;  if (n_bases < 2 * K + gap_min)
      return SYL_SUCCESS

   ;  acc_reset(cache)

                  /* initialize both dyad maps. Leave one base to be shifted in. */
   ;  for (l=0;l<K-1;l++)
      {  b1 = basemap[bases[l]]
      ;  hid1 |= b1 << 2*(K-2-l)
      ;  mask1 = mask1 << 1
      ;  if (b1 >= 4)
            mask1 |= 1
         ,  n_mask++

      ;  b2 = basemap[bases[l+K+gap_min]]
      ;  hid2 |= b2 << 2*(K-2-l)
      ;  mask2 = mask2 << 1
      ;  if (b2 >= 4)
         mask2 |= 1
   ;  }

                  /* if we print this, a spurious A fronts each component;
                   * it's ready to be shifted out
                  */
;if(0)fprintf(stderr, "start left %s right %s %d\n", get_sylmer(K, hid1, buf1), get_sylmer(K, hid2, buf2), (int) hid2);

      for (l=n_bases-K-gap_min;l<n_bases;l++)
      n_mask += basemap[bases[l]] >= 4 ? 1 : 0

                  /* garantuees space for at least one right dyad */
   ;  for (l=K-1;l < n_bases - K -gap_min;l++)
      {  unsigned b3, hid3, mask3, r3
      ;  b1 = basemap[bases[l]]
      ;  mask1 = mask1 << 1
      ;  hid1 = ((hid1 << 2) | b1) & LOMEGA(K)
      ;  if (b1 >= 4)
            mask1 |= 1
         ,  n_mask++

      ;  r = l+K+gap_min

      ;  b2 = basemap[bases[r]]
      ;  mask2 = mask2 << 1
      ;  hid2 = ((hid2 << 2) | b2) & LOMEGA(K)
      ;  if (b2 >= 4)
         mask2 |= 1

      ;  if (mask1 & MERMASK(K))
         {  n_unsyl += gap_max - gap_min + 1
;if(0)fprintf(stderr, "%s\n%.*s%s%.*s%s (skip 1)\n---\n"
,  bases
,  (int) l+1-K
,  space
,  get_sylmer(K, hid1, buf1)
,  (int) (r-l-K-1)
,  space
,  get_sylmer(K, hid2, buf2)
)
         ;  continue
      ;  }

         hid3 = hid2, mask3 = mask2, r3 = r

      ;  do
         {  unsigned location = (l + r3) / 2
         ;  r3++
;if(0)fprintf(stderr, "---\n%s\n%.*s%s%.*s%s%s\n"
,  bases
,  (int) l+1-K
,  space
,  get_sylmer(K, hid1, buf1)
,  (int) (r3-l-K-1)
,  space
,  get_sylmer(K, hid3, buf2)
,  mask3 & MERMASK(K) ? " (skip 2)" : ""
)
         ;  if (mask3 & MERMASK(K))
            n_unsyl++
         ;  else
            {  unsigned p, hidhid = (hid1 << (2*K)) | hid3
            ;  if (dst_arry)
               {  p = cache->acount->cells[hidhid].num
/* ,fprintf(stderr, "hidhid %s previous %d - %d\n", get_sylmer(KK, hidhid, buf1), (int) p, (int) location) */
               ;  if (!p || p+interval <= location || location+interval <= p)
                     dst_arry[hidhid].num++
                  ,  cache->acount->cells[hidhid].num = location
                  ,  n_syl++
               ;  else
                  n_r2++
/* ,fprintf(stderr, "skipping regular hidhid %s previous %d - %d\n", get_sylmer(KK, hidhid, buf1), (int) p, (int) location)
*/

            ;  }
               else
               {  unsigned o = -1u
               ;  HOFFSET_IF_PRESENT(o, cache->hcount, hidhid)
               ;  p = cache->hcount->cells[o].pl.num
/* ,fprintf(stderr, "hidhid %s location %d - %d\n", get_sylmer(KK, hidhid, buf1), (int) p, (int) location) */
               ;  if (o != -1u && (!p || p+interval <= location || location+interval <= p))
                     dst->hcount->cells[o].pl.num++
                  ,  cache->hcount->cells[o].pl.num = location
                  ,  n_syl++
               ;  else
                  n_r2++
/* ,fprintf(stderr, "skipping hash hidhid %s location %d - %d\n", get_sylmer(KK, hidhid, buf1), (int) p, (int) location)
*/

            ;  }
            }

            b3 = basemap[bases[r3]]
         ;  mask3 = mask3 << 1
         ;  hid3 = ((hid3 << 2) | b3) & LOMEGA(K)
         ;  if (b3 >= 4)
            mask3 |= 1
      ;  }
         while (r3 <= l + K+gap_max && r3 < n_bases)
;if(0)fprintf(stderr, "\n\n")
   ;  }
      dst->n_unit++
   ;  dst->n_bases +=   n_bases
   ;  dst->n_syl   +=   n_syl
   ;  dst->n_unsyl +=   n_unsyl
   ;  dst->n_r2    +=   n_r2
   ;  dst->n_mask  +=   n_mask
   ;  return SYL_SUCCESS
;  }


