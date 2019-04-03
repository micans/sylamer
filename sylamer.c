
   /* _____ Copyright 2008, Wellcome Trust Sanger Institute  ____________ */
   /* _____ Copyright 2008, 2009 European Bioinformatics Institute ______ */
   /* _____ Copyright 2010, 2011 European Bioinformatics Institute ______ */
   /* _____ Authors: Stijn van Dongen and Cei Abreu-Goodger  ____________ */
   /* _____ License: Gnu General Public License v3 or later  ____________ */
   /* ___________________________________________________________________ */


/* Programmer's contract.
 *    This code should be fully reentrant. It only employs (very) few constant
 *    global variables, and should free all memory in use before returning,
 *    both when halting prematurely (due to some error) and upon succesful
 *    completion.  Thus, sylmain() should cleanly load and execute from a
 *    shared object file.

 * TODO
 *    Further digest main(). Split off sibling programs.
 *
*/

/* Implementation note.
 *
 * .  Read query file, annotate with rank.
 * .  Sort queries so that they can be bsearched.
 * .  Read background fasta file, for each entry search queries.
 *       annotate fasta entry with rank in query if present, 
 *       otherwise with counter initialized at #query (subset mode).
 * .  Sort fasta entries based on annotated rank.
 *       this sorts the entries in query order.
 *
 * .  go through the sequence data twice.
 *       First time compute background.
 *       Second time compute windows/samples and optionally Markov conditionals.
 *
 * .  there is also code
 *    +  to only output the maximum over all bins,
 *    +  do repeated trials with random shuffling of utrs
 *    +  change count modes (-theshold and -unit-size)
 *    +  use Wilcox test rather than hyperg.
 *    +  etc etc (see DONE above and --help output).
 *
 * .  MAIN_WSTAT is a nice small example of how the interfaces work.
 * .  MAIN_COUNT as well, albeit not as clean.
 *
 *                General technique
 *
 * Encode k-mers in 2k bits, as unsigned integers.
 *
 * By default compute for all words.
 *    Then use perfect minimal hash by standard 4-ary alphabet map.  This
 *    perfect minimal hash simply acts as an offset into a fully utilized array
 *    (as all possible words are tracked).  Compute the hash efficiently, by
 *    utilizing the hash of the previous word and using bit shifting and
 *    masking.

 * With the -words and -read-expected options,
 *    use a real hash -- use an array that is smaller than the total possible
 *    number of words LSIZE(K). The idea is that we probably get more CPU cache
 *    hits that way for longer words -- as we assume that the word list is much
 *    smaller than LSIZE(K), the total number of possible words of that length.
 *    For words of size 9 for example, the full (integer) array of counts will
 *    have size 4 * 4**9 =~ 1M, exceeding usual cache sizes. Using a hash may
 *    thus significantly speed up the counting.  For smaller words, any speed
 *    losses measured have been minor.  Currently there is no way *not* to use
 *    a hash with the -words option.  However, we use the simplest possible way
 *    of hashing by taking the rightmost bits. So if the words file contains
 *    all possible words our hash values will be 0..LSIZE(K)-1 and the hash
 *    will be minimal and perfect so we hardly loose.

 CRUCIALNOTE
 *    The code works by generating a first hash, then cloning that hash
 *    for all other applications. This is crucial as we sometimes utilize
 *    the fact that each word is found at the same offset in every hash.

 *    The word list is kept, and is used to iterate through the hash arrays.
 *    This used to maintain word order as in the file, but these days
 *    we uniq and sort the words before hashing them.

 *                Input routines
 *
 * The line-reading routines are a bit messy/complicated. Some reasons:
 *    -  the low level read routine buffers using fread().
 *    -  support unix, mac, and dos line endings.
 *    -  support reading from gzipped files.
 *
 * Currently, the exact same code does input for unix, dos, and mac.  This
 * enables using any file on any system.  With that requirement (using any file
 * on any system), we have currently lost the ability to reliably detect
 * consecutive newlines (aka paragraph skips).  The DOS 'splice' event matters
 * here (introducing a paragraph skip where there was none), but there is
 * probably more nastiness lurking.  If there is ever a need to do this, the
 * code should #ifdef'd for different line-ending conventions, and the ability
 * to use files cross-platform should be dropped.

 *                Conventions
 * Some structure members and possibly variables have _p or _g suffix.
 * This is a leftover from the monolithic days, where _p indicated (user-adjustable)
 * parameter and _g indicated a global setting, state, or data object.
 * Some of these have been kept for silly reasons. For example, K_g and M_g
 * would just be so /short/ if it were just K and M.
*/


#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include "gsl/gsl_cdf.h"

#include "version.h"

#include "inttypes.h"
#include "interface.h"
#include "util.h"
#include "table.h"
#include "sylio.h"
#include "query.h"
#include "words.h"
#include "unit.h"
#include "array.h"
#include "acc.h"
#include "acciter.h"
#include "interval.h"


#ifndef USE_FLOAT
#define USE_FLOAT 0
#endif

#define FNAME_MAX_LENGTH 256


int exit_g = 0;


int mainfail
(  struct dataframe* dt
,  const char* fmt
,  ...
)
   {  if (fmt)
      {  va_list  args
      ;  va_start(args, fmt)
      ;  vfprintf(stderr, fmt, args)
      ;  fputc('\n', stderr)
      ;  va_end(args)
   ;  }
      dataframe_release(dt)
   ;  return exit_g ? exit_g : 1
;  }


unsigned get_mask
(  unsigned num
,  int K
,  dim hash_commit_fac
,  dim hash_commit_min
)
   {  unsigned mask

   ;  num = (hash_commit_fac * num) | hash_commit_min
   ;  mask = num

   ;  while ((num = num >> 1)) mask |= num
   ;  if (mask > LOMEGA(K))
      mask = LOMEGA(K)

   ;  return mask
;  }



   /* From a hash value that encodes a certain K-mer, compute the hash value
    * from an N-mer contained within that K-mer.
    *
    * The return result is some consecutive set of bits from sylid. Note that
    * we have to right-shift those bits so that the right-most bit is
    * encoded at the right most position of a C unsigned.
    *
    *     |__length_| (l)
    * + + + + + + + + + + +  K
    *     ^offset (o)
    * 0 1 2 3 4 5 6 7 8 9 10
    *
    * We have to shift K-l-o bits. As we do things 4-ary, it is multiplied by two.
    *
   */
unsigned get_submer_id
(  unsigned sylid          /* encodes k-mer of length K           */
,  unsigned K              /* length K                            */
,  unsigned offset         /* start at offset                     */
,  unsigned m              /* the markov word size we are testing */
)
   {  return ((sylid >> 2 * (K - offset - m)) & LOMEGA(m))
;  }


void dump_expect
(  FILE* fp
,  const char* word
,  unsigned long i
,  double e
)
   {  double ey = get_string_entropy(word)
   ;  fprintf(fp, "%s\t%.8g\t%.4g\t%lu\n", word, e, ey, i)
;  }



sylstatus get_mkv_expected
(  unsigned hid
,  const rcell* SA            /* the smaller word size (nominator) */
,  const rcell* SB            /* the larger word size (numerator)  */
,  dim SA_total
,  dim SB_total
,  dim K
,  dim m                      /* m == 1 equivalent with SA_total == 0 */
,  double* valp
)
   {  double nom = SA_total ? 1.0 : 0.0
   ;  double num = SB_total ? 1.0 : 0.0
   ;  char buf[16]
   ;  dim j

;if(0)fprintf(stderr, "num %g\n", num)
   ;  if (num)
      for (j=0;j<=K-m;j++)
      {  unsigned ii = get_submer_id(hid, K, j, m)
      ;  num *= SB[ii].num * 1.0 / SB_total
;if(0)fprintf(stderr, "num %g\n", num)
   ;  }
      if (m>1 && nom)
      for (j=1;j<=K-m;j++)
      {  unsigned ii = get_submer_id(hid, K, j, m-1)
      ;  nom *= SA[ii].num * 1.0 / SA_total
   ;  }

      if (m == 1)
      {  if (!num)
         {  miaow
            (  "panic: no numerator at sylid %u/%s with numerator %g (should be impossible)"
            ,  hid
            ,  get_sylmer(K, hid, buf)
            ,  num
            )
         ;  return SYL_FAIL
      ;  }
         *valp = num
   ;  }
      else
      *valp = num/nom
   ;  return SYL_SUCCESS   
;  }


   /*
   */
sylstatus have_expectations_gapped
(  const char* prefix
,  dim   binid             /* 0 for background    */
,  dim   K
,  dim   m
,  dim   n_bases
,  dim   gap_size          /*                     */
,  const rcell* SA         /* size LSIZE(K/2) - 1 */
,  const rcell* SB         /* size LSIZE(K/2) - 0 */
,  mkv*  markov
,  int use_markov_as_usual
)
   {  char buf[16]
   ;  char fname[FNAME_MAX_LENGTH]
   ;  FILE* iffpout = NULL
   ;  sylstatus status = SYL_SUCCESS
   ;  ofs o = 0
   ;  unsigned hid = 0
   ;  float* ep = NULL
   ;  unsigned Kh = K/2
   ;  char buf1[16]
   ;  char buf2[16]
   ;  dim SA_total = 0
   ;  dim SB_total = 0
   ;  dim i

   ;  if (use_markov_as_usual)
      {  for (i=0;i<LSIZE(m-1);i++)
         SA_total += SA[i].num
      ;  for (i=0;i<LSIZE(m-0);i++)
         SB_total += SB[i].num
   ;  }

      if (prefix)
      {  if (binid)
         myfname(fname, FNAME_MAX_LENGTH, "%s.%u", prefix, (unsigned) binid)
      ;  else
         myfname(fname, FNAME_MAX_LENGTH, "%s.bg", prefix)
      ;  iffpout = myfopen(fname, "w")         /* subsequent use is conditional */
   ;  }

      while (0 <= (o = mkv_iterate(o, markov, &hid, &ep, &status)))
      {  unsigned hid1 = (hid >> K) & LOMEGA(Kh)
      ;  unsigned hid2 = hid & LOMEGA(Kh)

      ;  if (use_markov_as_usual)
         {  double val1 = 0.0, val2 = 0.0
         ;  if ((status = get_mkv_expected(hid1, SA, SB, SA_total, SB_total, Kh, m, &val1)))
            break
         ;  if ((status = get_mkv_expected(hid2, SA, SB, SA_total, SB_total, Kh, m, &val2)))
            break
         ;  *ep = val1 * val2 * gap_size
      ;  }
         else
         *ep = (SB[hid1].num * 1.0 * SB[hid2].num * 1.0 * gap_size) / (n_bases * 1.0 * n_bases)

;if(0)fprintf(stdout, "%2d->%s %2d->%s %d %d %d %g\n",
(int) hid1, get_sylmer(Kh, hid1, buf1),
(int) hid2, get_sylmer(Kh, hid2, buf2),
(int) n_bases,
(int) SB[hid1].num, (int) SB[hid2].num,
(double) (*ep * n_bases)
)
      ;  if (iffpout)
         dump_expect(iffpout, get_sylmer(K, hid, buf), hid, *ep)
   ;  }

      if (iffpout)
      myfclose(&iffpout)
   ;  return status
;  }




   /* Compute expected frequencies of occurrence for each word of length K,
    * based on the (normalized) frequencies of occurrence for words of length m
    * and m-1.  It is computed as the product of frequencies of all subwords of
    * length m, divided by the product of frequencies of all subwords of length
    * m-1 EXCEPT the leading and trailing word of length m-1.
   */
sylstatus have_expectations
(  const char* prefix
,  dim   binid        /* 0 for background */
,  dim   K
,  dim   m
,  rcell* SA          /* size LSIZE(m-1) */
,  rcell* SB          /* size LSIZE(m)   */
,  mkv*  markov
)
   {  dim SA_total = 0
   ;  dim SB_total = 0
   ;  char buf[16]
   ;  char fname[FNAME_MAX_LENGTH]
   ;  FILE* iffpout = NULL
   ;  sylstatus status = SYL_SUCCESS
   ;  ofs o = 0
   ;  unsigned hid = 0
   ;  float* ep = NULL
   ;  dim i

   ;  if (prefix)
      {  if (binid)
         myfname(fname, FNAME_MAX_LENGTH, "%s.%u", prefix, (unsigned) binid)
      ;  else
         myfname(fname, FNAME_MAX_LENGTH, "%s.bg", prefix)
      ;  iffpout = myfopen(fname, "w")         /* subsequent use is conditional */
   ;  }

      for (i=0;i<LSIZE(m-1);i++)
      SA_total += SA[i].num
   ;  for (i=0;i<LSIZE(m-0);i++)
      SB_total += SB[i].num            /* hv 2 SB_total == 0 */

#if 0
,fprintf(stderr, ". %d %u\n", (int) i, (int) SB[i].num)
#endif

   ;  while (0 <= (o = mkv_iterate(o, markov, &hid, &ep, &status)))
      {  double val
      ;  if ((status = get_mkv_expected(hid, SA, SB, SA_total, SB_total, K, m, &val)))
         break
      ;  *ep = val

;if(0)fprintf(stderr, "%u expect %f\n", (unsigned) m, (double) *ep)
      ;  if (iffpout)
         dump_expect(iffpout, get_sylmer(K, hid, buf), hid, *ep)
   ;  }

      if (iffpout)
      myfclose(&iffpout)
   ;  return status
;  }


int myfinish (struct dataframe* dt, int verbose)
   {
if (verbose_G && verbose) {
fprintf(stderr, "\nPlease cite:\n");
fprintf(stderr, "    Stijn van Dongen, Cei Abreu-Goodger & Anton J. Enright,\n");
fprintf(stderr, "    Detecting microRNA binding and siRNA off-target effects\n");
fprintf(stderr, "    from expression data, Nature Methods 2008 (PMID 18978784).\n");
fputc('\n', stderr);
}
      dataframe_release(dt)
   ;  return SYL_SUCCESS
;  }


#define INTERFACE(par, mode)  (par.modes & (mode))
#define INTERFACE2(par, mode) (par.modes2 & (mode))      /* fixme (programmer-unfriendly) */
#define MAIN_MODE(par, mode)  (par.main_modes & (mode))


void output_utrtable
(  struct parameter* pp
,  struct dataframe* df
,  dim n_chunks
)
   {  dim i, j
   ;  char sylbuf[16]
   ;  for (i=0;i<df->wordlist.n_words;i++)
      {  fputs(get_sylmer(pp->K_g, df->wordlist.words ? df->wordlist.words[i] : i, sylbuf), stdout)
      ;  for (j=i*n_chunks;j<(i+1)*n_chunks;j++)
         fprintf(stdout, "\t%g", df->utrtable[j])
      ;  fputc('\n', stdout)
   ;  }
   }


void output_table
(  struct parameter* pp
,  struct dataframe* df
,  dim trial_i
)
   {  dim i, j
   ;  char sylbuf[16]
   ;  dim n_bin = df->n_bin
   ;  dim n_bin_used = pp->lastbin && pp->lastbin < n_bin ? pp->lastbin : n_bin
   ;  dim n_peaktrough = 0
   ;  unsigned char* select_table
      =     pp->toptable
         ?  get_toptable(df->table, df->wordlist.n_words, df->n_bin, n_bin_used, pp->toptable)
                  /* get_toptable might return NULL, subsequent use conditional */
         :  NULL

   ;  if (!trial_i)
      {  fputs("upper", pp->fpresult)
      ;  for (i=0;i<n_bin_used;i++)
         fprintf(pp->fpresult, "\t%lu", (unsigned long) df->bins_info[i].offset)
      ;  fputc('\n', pp->fpresult)
      ;  if (INTERFACE2(pp[0], Mode_UTRATTR))
         {  fputs("utrlen", pp->fpresult)
         ;  for (i=0;i<n_bin_used;i++)
            fprintf(pp->fpresult, "\t%5g", (double) (df->bins_info[i].utr_length * 1.0 / df->bins_info[i].size))
         ;  fputc('\n', pp->fpresult)

         ;  fputs("nmasked", pp->fpresult)
         ;  for (i=0;i<n_bin_used;i++)
            fprintf(pp->fpresult, "\t%5g", (double) (df->bins_info[i].n_masked * 1.0 / df->bins_info[i].size))
         ;  fputc('\n', pp->fpresult)

         ;  fputs("gccontent", pp->fpresult)
         ;  for (i=0;i<n_bin_used;i++)
            fprintf(pp->fpresult, "\t%5g", (double) (df->bins_info[i].gccontent * 1.0 / df->bins_info[i].size))
         ;  fputc('\n', pp->fpresult)
      ;  }
      }

      for (i=0;i<df->wordlist.n_words;i++)
      {  if (select_table && !select_table[i])
         continue
      ;  if (select_table && select_table[i] == 3)
         n_peaktrough++
      ;  fputs(get_sylmer(pp->K_g, df->wordlist.words ? df->wordlist.words[i] : i, sylbuf), pp->fpresult)
      ;  if (pp->n_trial > 1)
         fprintf(pp->fpresult, ".%d", (int) trial_i)
      ;  for (j=i*n_bin;j<i*n_bin+n_bin_used;j++)
         fprintf(pp->fpresult, "\t%g", df->table[j])
      ;  fputc('\n', pp->fpresult)
   ;  }

      if (select_table)
      {  tell(VERBOSE_ONCE, "%lu elements intersect topmin topmax", n_peaktrough)
      ;  free(select_table)
   ;  }
   }


typedef struct rank_unit
{  long     index
;  double   ord
;  double   value
;
}  rank_unit   ;


void*  rank_unit_init
(  void*   ruv
)
   {  rank_unit* ru =   (rank_unit*) ruv
   ;  ru->index     =   -1
   ;  ru->value     =   -1
   ;  ru->ord       =   -1.0
   ;  return ru
;  }


int rank_unit_cmp_index
(  const void* rua
,  const void* rub
)
   {  long a = ((rank_unit*) rua)->index
   ;  long b = ((rank_unit*) rub)->index
   ;  return a < b ? -1 : a > b ? 1 : 0
;  }


int rank_unit_cmp_value
(  const void* rua
,  const void* rub
)
   {  double a = ((rank_unit*) rua)->value
   ;  double b = ((rank_unit*) rub)->value
   ;  return a < b ? 1 : a > b ? -1 : 0
;  }


   /* This uses df->utrtable, which contains per-word/per-utr results,
    * in word-major order.
    * One of the main differences with hypergeometric code is that here
    * we loop over the words first, rather than the binned UTRs.
   */
static sylstatus compute_mww
(  struct parameter* pp
,  struct dataframe* df
)
   {  sylstatus status = SYL_SUCCESS
   ;  dim n_analyse = df->n_analyse
   ;  dim n_bin = df->n_bin
   ;  rank_unit* rus = mycalloc((1+n_analyse), sizeof rus[0])
   ;  dim w

   ;  if (!rus)
      return arrr("could not alloc rus array")

   ;  for (w=0;w<df->wordlist.n_words;w++)
      {  dim i, j, stretch = 0
      ;  for (j=0;j<n_analyse;j++)
         {  rus[j].index = j
         ;  rus[j].value = df->utrtable[w*n_analyse+j]
         ;  rus[j].ord = -1.0
      ;  }

         qsort(rus, n_analyse, sizeof rus[0], rank_unit_cmp_value)

      ;  rus[n_analyse].value = rus[n_analyse-1].value + 1000  /* sentinel out-of-bounds value */
      ;  stretch = 1
      ;  for (i=0;i<n_analyse;i++)
         {  if (rus[i].value != rus[i+1].value)           /* input current stretch */
            {   for (j=0;j<stretch;j++)
                rus[i-j].ord = 1 + (i+1 + i-stretch) / 2.0
            ;   stretch = 1
        ;   }
            else
            stretch++
      ;  }
         qsort(rus, n_analyse, sizeof rus[0], rank_unit_cmp_index)

;if(0)
for(i=0;i<n_analyse;i++)
fprintf(stderr, "word %d value %g idx %d\n", (int) w, rus[i].value, (int) rus[i].ord)

      ;  {  long n1 = 0, n2 = n_analyse
         ;  double R1 = 0.0
         ;  dim x
         ;  spacer sp = { 0 }             /* initialisation required right here */

         ;  for (x=0;x<n_bin;x++)
            {  double sign, U, sigma, mean, p
            ;  if (pp->lastbin && x >= pp->lastbin)
               break
            ;  interval_step(df->thecollection, df->n_analyse, &sp)

            ;  sp.i_hi_prev = sp.i_hi

            ;  if (!w)
               tell
               (  VERBOSE_BIN
               ,  "bin %lu [0..%lu) %lu"
               ,  (unsigned long) (x+1)
               ,  (unsigned long) sp.i_hi
               )

            ;  if (!w && df->bins_info)
               {  df->bins_info[x].offset    =  sp.i_hi
               ;  df->bins_info[x].size      =  sp.i_hi - sp.i_lo
               ;  df->bins_info[x].utr_length = 0.0
               ;  df->bins_info[x].n_masked  = 0.0
               ;  df->bins_info[x].gccontent = 0.0
            ;  }

               for (i=sp.i_lo;i<sp.i_hi;i++)
               {  n1++; n2--
               ;  R1 += rus[i].ord
               ;  if (!w && df->bins_info)
                  {  df->bins_info[x].utr_length += df->thecollection->list[i].n_bases
                  ;  df->bins_info[x].n_masked   += df->thecollection->list[i].n_masked
                  ;  df->bins_info[x].gccontent  += df->thecollection->list[i].gccontent
               ;  }
               }
               mean = n1 * 0.5 * n2
            ;  sigma = sqrt(n1 * 1.0 * n2 * 1.0 * (n1 + n2 + 1.0) / 12.0)
            ;  U = R1 - (n1 * (n1 + 1.0)) * 0.5
            ;  sign = U < mean ? -1 : 1
;if(0)fprintf(stderr, "%g %g %g %g", mean, U, rus[i].value, rus[i].ord)

            ;  if (sign && sigma)
               p =   U < mean
                  ?  gsl_cdf_ugaussian_P(( U - mean ) / sigma)
                  :  gsl_cdf_ugaussian_Q(( U - mean ) / sigma)
            ;  else
               p = 1.0

            ;  df->table[w * n_bin + x] = (p ? log(p) / log(10) : 312) * sign
;if(0)fprintf(stderr, " %g\n", df->table[w * n_bin + x])
         ;  }

            if (!pp->lastbin && (n2 || n1 != n_analyse))
            {  fprintf(stderr, "mww-U fail n2=%d n1=%d\n", (int) n2, (int) n1)
            ;  status = SYL_FAIL
            ;  goto DONE
         ;  }
         }
      }

   DONE:
      free(rus)
   ;  return status
;  }


            /* dt.table is malloced by caller.
             * Do we need all the acc_init calls prior to our invocation?
             * use pp->binsize
             * TODO: separate more cleanly the case word_stat
             *    tidy up init, ownership, dimensions.
            */
static sylstatus dispatch_mww
(  struct parameter* pp
,  struct dataframe* df
,  int word_stat
)
   {  acc ctwindow
   ;  char sylbuf[16]
   ;  mkv* markov_b = pp->M_g ? &df->mkv1 : NULL
   ;  mkv* markov_l = pp->M_g ? &df->mkv2 : NULL
   ;  dim n_analyse = df->n_analyse              /* fixme docme why not n_analyse */
   ;  unsigned r2check = INTERFACE(pp[0], MODE_R2ON)
   ;  sylstatus status = SYL_SUCCESS
   ;  acc_iter iterator
   ;  dim x

   ;  if (!word_stat)
      {  dim ncalloc =  n_analyse * df->wordlist.n_words
      ;  if (!(df->utrtable = mycalloc(ncalloc, sizeof df->utrtable[0])))
         return arrr("could not alloc utr table")
;fprintf(stderr, "ncalloc %d, analyse %d universe %d lastbin %d\n", (int) ncalloc, (int) n_analyse, (int) df->n_universe, (int) pp->lastbin)
   ;  }

      if (acc_init(&ctwindow, pp->K_g, df->hash_template_g))
      return arrr("acc_init failure")

   ;  if (pp->M_g)
      {  dim x
      ;  for (x=0;x<n_analyse;x++)
         {  unit* ut = df->thecollection->list+x
         ;  if (acc_addto(pp->M_g,   ut, &df->bg_xyz, r2check))
            return arrr("addto xyz background")
         ;  if (acc_addto(pp->M_g-1, ut, &df->bg_xy , r2check))
            return arrr("addto xy background")
      ;  }
         if
         (  have_expectations
            (  NULL
            ,  0                 /* signifies background */
            ,  pp->K_g
            ,  pp->M_g
            ,  df->bg_xy.acount->cells          /* tempdebug these seem to be zero */
            ,  df->bg_xyz.acount->cells         /* tempdebug these seem to be zero */
            ,  markov_b
         )  )
         return arrr("failed to compute background expectations")
   ;  }


      if (word_stat)
      {  int i
      ;  fputs("id\tnbases", pp->fpresult)
      ;  for (i=0;i<df->wordlist.n_words;i++)
            fputc('\t', pp->fpresult)
         ,  fputs(get_sylmer(pp->K_g, df->wordlist.words ? df->wordlist.words[i] : i, sylbuf), pp->fpresult)
      ;  fputc('\n', pp->fpresult)
   ;  }

      for (x=0;x<n_analyse;x++)
      {  unit* ut = df->thecollection->list+x
      ;  if (word_stat && pp->lastbin && x >= pp->lastbin)
         break
                     /* ^ we need to compute everything in order to sort and resort
                      * later (in the case of mww), so cannot skip computation here.
                     */

      ;  acc_reset(&ctwindow)
      ;  if (acc_addto(pp->K_g, ut, &ctwindow, MODE_R2ON))
         return arrr("addto window")

      ;  if (word_stat)
         fprintf(pp->fpresult, "%s\t%lu", ut->sid, (ulong) ut->n_bases)

      ;  if (pp->M_g)
         {  acc_reset(&df->window_xy)
         ;  acc_reset(&df->window_xyz)
         ;  if (acc_addto(pp->M_g,  ut, &df->window_xyz, MODE_R2ON))
            return arrr("addto xyz window")
         ;  if (acc_addto(pp->M_g-1, ut, &df->window_xy , MODE_R2ON))
            return arrr("addto xy window")
         ;  if
            (  have_expectations
               (  NULL
               ,  1                 /* signifies not-background, ignored otherwise */
               ,  pp->K_g
               ,  pp->M_g
               ,  df->window_xy.acount->cells
               ,  df->window_xyz.acount->cells
               ,  markov_l
            )  )
            return arrr("failed to compute window expectations")
      ;  }

         acc_iter_reset(&iterator)

      ;  while (acc_bi_iterate(&iterator, &df->bg, &ctwindow, markov_b, markov_l, df->eoracle_g, &status))
         {  dim window_count_y=  iterator.count2
         ;  float elocal_y    =  iterator.val2
         ;  dim y =  iterator.hid
         ;  double freq = ctwindow.n_syl ? (window_count_y * 1.0) / ctwindow.n_syl : 0.0
         ;  double ratio = pp->M_g && elocal_y ? freq / elocal_y : 1
         ;  double result
            =  pp->M_g
               ?  (  INTERFACE(pp[0], MODE_LOGFOLD)
                  ?  ( freq ? log(ratio) / log(10.0) : -312 )
                  :  freq - elocal_y
                  )
               :  freq
         ;  if (!result && INTERFACE(pp[0], MODE_PERTURB_TIES))
            result = ((random() >> 10) * 0.001) / RAND_MAX
                  /*
                   * fixme considerme: for freq - local_y
                   * one might have a range of [-delta, delta]
                   * (so, including negative numbers). At
                   * the moment it is [0, delta], where delta
                   * is smaller than the smallest possible
                   * frequency with high probability - exceptions
                   * require large UTR length and lucky roll of die.
                  */

         ;  if (word_stat)
            fprintf(pp->fpresult, "\t%.5g", result)
         ;  else if (1)
            df->utrtable[(iterator.offset-1) * n_analyse + x] = result
;if(0)fprintf
   (  pp->fpresult
   ,  "%s\t%lu\t%5g\t%5g\t%lu\t%lu\n"
   ,  get_sylmer(pp->K_g, y, sylbuf)
   ,  y
   ,  elocal_y
   ,  ctwindow.n_syl ? (window_count_y * 1.0) / ctwindow.n_syl : 0
   ,  window_count_y
   ,  ctwindow.n_syl
   )
      ;  }

         if (word_stat)
         fputc('\n', pp->fpresult)

      ;  if (status)
         return arrr("mww error while iterating over counts")
   ;  }

      if (df->utrtable)
      {  if (0)
         output_utrtable(pp, df, n_analyse)
      ;  if (compute_mww(pp, df))
         return arrr("mww error while computing mww")
      ;  output_table(pp, df, 0)
   ;  }

            /* fixme; these free/close instances are not yet designed to follow error paths
             * (i.e. mainfail()).
             * Solution: put them together in a new data bundle and introduce
             * a callback routine in dataframe struct.
            */
      myfclose(&pp->fpresult)
   ;  acc_free(&ctwindow)
;fprintf(stderr, "hero\n")
   ;  return SYL_SUCCESS
;  }



static sylstatus dispatch_count
(  struct parameter* pp
,  struct dataframe* df
)
   {  acc ctwindow, ctthreshold, ctcache
   ;  acc* ctacc_threshold = pp->per_sequence || pp->cap ? &ctthreshold : NULL
   ;  dim n_analyse = df->n_analyse
   ;  sylstatus status = SYL_SUCCESS
   ;  acc_iter iterator
   ;  dim tblidx = 0, i
                                                /* 4 extra because of ctwindow.n_***
                                                 * hack further below.
                                                */
   ;  if (!(df->tablect = mycalloc(n_analyse * (df->wordlist.n_words+4), sizeof df->tablect[0])))
      return mainfail(df, "could not alloc table")

   ;  if (acc_init(&ctwindow, pp->K_g, df->hash_template_g))
      return mainfail(df, NULL)

   ;  if (INTERFACE(pp[0], MODE_GAPPY) && acc_init(&ctcache, pp->K_g, df->hash_template_g))
      return mainfail(df, NULL)

   ;  if (ctacc_threshold && acc_init(ctacc_threshold, pp->K_g, df->hash_template_g))
      return mainfail(df, NULL)

   ;  for (i=0;i<n_analyse;i++)
      {  unit* ut = df->thecollection->list+i
      ;  if (pp->lastbin && i >= pp->lastbin)
         break

      ;  acc_reset(&ctwindow)

      ;  status
         =  INTERFACE(pp[0], MODE_GAPPY)
         ?  acc_addto_gapped(pp->K_g, ut, &ctwindow, pp->gap_min, pp->gap_max, pp->gap_gap, &ctcache)
         :  acc_addto(pp->K_g, ut, &ctwindow, INTERFACE(pp[0], MODE_R2ON))
      ;  if (status)
         return mainfail(df, "addto count")

      ;  acc_iter_reset(&iterator)

      ;  while (acc_iterate(&iterator, &ctwindow, &status)) 
         df->tablect[tblidx++]  =  iterator.count1

      ;  if (status)
         return mainfail(df, "error while iterating over counts")

      ;  df->tablect[tblidx++]   =  ctwindow.n_bases
      ;  df->tablect[tblidx++]   =  ctwindow.n_mask
      ;  df->tablect[tblidx++]   =  ctwindow.n_syl
      ;  df->tablect[tblidx++]   =  ctwindow.n_unsyl
   ;  }

      {  char sylbuf[16]
      ;  dim n_bin_used = pp->lastbin && pp->lastbin < n_analyse ? pp->lastbin : n_analyse
      ;  dim n_passed = 0
      ;  if (MAIN_MODE(pp[0], MAIN_FILTER))
         {  for (i=0;i<n_bin_used;i++)
            {  dim j, sum = 0
            ;  int ok = 1
            ;  for (j=0;j<df->wordlist.n_words;j++)
               {  dim ct = df->tablect[i*(df->wordlist.n_words+4) + j]
               ;  if (INTERFACE(pp[0], MODE_SUM))
                  sum += ct
               ;  else
                  {  ok =
                     (  (pp->unit_at_most  < 0 || ct <= pp->unit_at_most )
                     && (pp->unit_at_least < 0 || ct >= pp->unit_at_least)
                     )  ?
                     1  :  0
                  ;  if ( (!ok && INTERFACE(pp[0], MODE_AND)) || (ok && !INTERFACE(pp[0], MODE_AND)) )
                     break
               ;  }
               }
               if (INTERFACE(pp[0], MODE_SUM))
               ok =
               (  (pp->unit_at_most  < 0 || sum <= pp->unit_at_most )
               && (pp->unit_at_least < 0 || sum >= pp->unit_at_least)
               )  ?
               1  :  0
            ;  if (ok)
                  fputs(df->thecollection->list[i].sid, pp->fpresult)
               ,  fprintf(pp->fpresult, "\t%u\n", (unsigned) (i+1))
               ,  n_passed++
         ;  }
            tell(VERBOSE_ONCE, "sequences passed %lu", (unsigned long) n_passed)
      ;  }
         else if (MAIN_MODE(pp[0], MAIN_COUNT))
         {  fputs("id", pp->fpresult)
         ;  for (i=0;i<df->wordlist.n_words;i++)
               fputc('\t', pp->fpresult)
            ,  fputs(get_sylmer(pp->K_g, df->wordlist.words ? df->wordlist.words[i] : i, sylbuf), pp->fpresult)
         ;  fputs("\tbases\tmasked\tsyl\tunsyl\trank\n", pp->fpresult)
         ;  for (i=0;i<n_bin_used;i++)
            {  dim j
            ;  fputs(df->thecollection->list[i].sid, pp->fpresult)
            ;  for (j=0;j<df->wordlist.n_words+4;j++)
               fprintf(pp->fpresult, "\t%lu", df->tablect[i*(df->wordlist.n_words+4) + j])
            ;  fprintf(pp->fpresult, "\t%lu", (unsigned long) (i+1))
            ;  fputc('\n', pp->fpresult)
         ;  }
         }
      }
            /* fixme; these free/close instances are not yet designed to follow error paths
             * solution: put them together in a new data bundle and introduce
             * a callback routine in dataframe struct.
            */
      myfclose(&pp->fpresult)
   ;  acc_free(&ctwindow)
   ;  if (INTERFACE(pp[0], MODE_GAPPY))
      acc_free(&ctcache)
   ;  if (ctacc_threshold)
      acc_free(ctacc_threshold)
   ;  return myfinish(df, 0)
;  }


static sylstatus dispatch_wstat        /* mq */
(  struct parameter* pp
,  struct dataframe* df
)
   {  char sylbuf[16]
   ;  sylstatus status = SYL_SUCCESS
   ;  acc_iter iterator
   ;  mkv* markov_b = pp->M_g ? &df->mkv1 : NULL
   ;  mkv* markov_l = pp->M_g ? &df->mkv2 : NULL
   ;  dim bg_count  = df->bg_use_count
   ;  dim j

   ;  fprintf(pp->fpresult, "word\twordid\tfreqbg\tfreqxp\tfreqxpa\n")

   ;  if (!pp->M_g)
      return mainfail(df, "-m parameter required for --word-stat")

                              /* hv 1 check this (added to fix hv 2) */
                              /* fixme docme why is this window rather than background? */
   ;  for (j=0;j<df->n_analyse;j++)
      {  unit* ut = df->thecollection->list+j
      ;  if (acc_addto(pp->K_g, ut, &df->window, INTERFACE(pp[0], MODE_R2ON)))
         return mainfail(df, "addto window")
      ;  if (acc_addto(pp->M_g,  ut, &df->window_xyz, INTERFACE(pp[0], MODE_R2ON)))
         return mainfail(df, "addto xyz window")
      ;  if (acc_addto(pp->M_g-1, ut, &df->window_xy , INTERFACE(pp[0], MODE_R2ON)))
         return mainfail(df, "addto xy window")
   ;  }

      if                      /* fixme: document this and below have_expectation invocation */
      (  have_expectations
         (  NULL
         ,  1                 /* signifies not-background, ignored otherwise */
         ,  pp->K_g
         ,  pp->M_g
         ,  df->window_xy.acount->cells
         ,  df->window_xyz.acount->cells
         ,  markov_l
      )  )
      return mainfail(df, "failed to compute window expectations")

   ;  if
      (  have_expectations
         (  NULL
         ,  0                 /* signifies background */
         ,  pp->K_g
         ,  pp->M_g
         ,  df->bg_xy.acount->cells
         ,  df->bg_xyz.acount->cells
         ,  markov_b
      )  )
      return mainfail(df, "failed to compute background expectations")

   ;  acc_iter_reset(&iterator)

   ;  while (acc_bi_iterate(&iterator, &df->bg, &df->window, markov_b, markov_l, 0, &status))
      {  dim bg_count_y    =  iterator.count1
      ;  float eglobal_y   =  iterator.val1
      ;  float elocal_y    =  iterator.val2

      ;  double anchor_fac =  bg_count_y * 1.0 / (1.0 * eglobal_y * bg_count)
      ;  dim y =  iterator.hid            /* should never be used as offset now */

      ;  fprintf
         (  pp->fpresult
         ,  "%s\t%lu\t%.4g\t%.4g\t%.4g\n"
         ,  get_sylmer(pp->K_g, y, sylbuf)
         ,  y
         ,  (double) bg_count_y * 1.0 / bg_count
         ,  (double) elocal_y
         ,  (double) elocal_y * anchor_fac
         )
   ;  }
      myfclose(&(pp->fpresult))
   ;  return myfinish(df, status)
;  }


int sylmain
(  int argc
,  char* argv[]
)

   {  query             qry
   ;  collection        col
   ;  struct parameter  par            /* mostly encodes user settings and derived parameters */
   ;  struct dataframe  dt             /* encodes the bundle of data that we work with */

   ;  dim trial_i = 0
   ;  dim i = 0

   ;  time_t t = time(NULL)
   ;  pid_t p  = getpid()
   ;  pid_t pp = getppid()
   ;  sylstatus status = SYL_SUCCESS

   ;  unsigned long s
      =     (p ^ p << 4 ^ p << 16 ^ p << 28)
         ^  (pp ^ pp << 8 ^ pp << 24)
         ^  (t ^ t << 12 ^ t << 20)
   ;  srandom(s)

   ;  mod_acc_init()
   ;  query_init(&qry)
   ;  collection_init(&col)
   ;  parameter_init(&par)
   ;  dataframe_null(&dt, &par, &col, &qry)

   ;  status = parse_args(&par, argc, argv)

                                       /* fixme (trivially) unify myfinish/mainfail */
   ;  if (status == SYL_DONE)
      return myfinish(&dt, 0)
   ;  else if (status == SYL_FAIL)
      return mainfail(&dt, NULL)


   ;  if (par.fnwords && !(par.fpwords = myfopen(par.fnwords, "r")))
      return mainfail(&dt, NULL)

   ;  if (!par.fpfasta)
      return mainfail(&dt, "need fasta input file (see -h for options)")

   ;  if (par.patlist && INTERFACE(par, MODE_FREQUENCIES))
      return mainfail(&dt, "-read-expected and -w do not rumba")


/* ***************************
  *   move this to words.c, probably
   */

   ;  if (par.fpwords || par.patlist || par.xpatlist || par.dyad_modes || par.entropy)
      {  do
         {  if (par.dyad_modes)
            {  if (!par.K_g)
               par.K_g = 6
            ;  par.K_g |= 1, par.K_g ^= 1
            ;  status = get_dyadlist(par.dyad_modes, &dt.wordlist, par.K_g >> 1)
         ;  }
            else if (par.fpwords)
            {  status
               =  get_wordlist
                  (  par.fpwords
                  ,  &dt.wordlist
                  ,  INTERFACE(par, MODE_FREQUENCIES) ? &dt.eoracle_g : NULL
                  ,  &par.K_g
                  )
            ;  myfzclose(&par.fpwords)
         ;  }
            if (status)
            return mainfail(&dt, "error building words from %s", par.dyad_modes ? "dyads" : "file")

         ;  if (INTERFACE(par, MODE_FREQUENCIES))
            {  if (par.patlist || par.xpatlist || par.entropy)
               miaow("options -w/-xw/-entropy-gq do not work with -read-expected")
            ;  break
         ;  }

            if
            (  par.patlist
            && (status = get_bypattern(par.patlist, &dt.wordlist, &par.K_g))
            )
            return mainfail(&dt, "error building -w supplied words")

         ;  if (!dt.wordlist.words && (status = wlist_full(&dt.wordlist, par.K_g)))
            return mainfail(&dt, "error building full word list")

         ;  if (par.xpatlist)
            {  wlist wdelete = { NULL, 0, 0 }
            ;  if ((status = get_bypattern(par.xpatlist, &wdelete, &par.K_g)))
               return mainfail(&dt, "error building -xw supplied words")
            ;  wlist_remove(&dt.wordlist, &wdelete)
            ;  wlist_free(&wdelete)
         ;  }

            if (par.entropy)
            {  dim n_list = dt.wordlist.n_words
            ;  dim n_entr = wlist_entropy(&dt.wordlist, par.K_g, par.entropy)
            ;  tell
               (  VERBOSE_ONCE
               ,  "removed %d from word list of size %d by entropy filter %.6f"
               ,  (int) n_entr
               ,  (int) n_list
               ,  par.entropy
               )
         ;  }

            if (par.modes & MODE_BAG)
            get_bags(&dt.wordlist, par.K_g)
      ;  }
         while (0)
   ;  }
      else
      {  if (!par.K_g)
         par.K_g = 6
      ;  dt.wordlist.n_words = LSIZE(par.K_g)
               /* fixme docme: why is it not necessary to touch dt.wordlist.words?
                * (quite likely it's NULL, but n_words is overloaded .. sigh)
               */
   ;  }


   /* ***************************
  *   move above to words.c, probably
 */

               /* this is both dubious and stress-test-useful: we use dt.wordlist.words as the
                * test whether to use the wordlist, even if it is empty (that is,
                * dt.wordlist.n_words == 0, e.g. by entropy removal or word exclusion).
               */
      if (dt.wordlist.words)
      {  if (!dt.wordlist.n_words)
         miaow("no words! that is a unit test, that is, hold steady")
      ;  else if (dt.wordlist.n_words)
         {  dim n_list = dt.wordlist.n_words
         ;  dim n_dup = wlist_uniq(&dt.wordlist)
         ;  if (n_dup)
            tell
            (  VERBOSE_ONCE
            ,  "removed %d duplicates from word list of size %d"
            ,  (int) n_dup
            ,  (int) n_list
            )
      ;  }
         dt.hmask_g = get_mask(dt.wordlist.n_words, par.K_g, par.hash_commit_fac_g, par.hash_commit_min_g)
      ;  dt.hash_template_g = harray_new(&dt.wordlist,  dt.hmask_g, par.K_g)
      ;  if (!dt.hash_template_g)
         return mainfail(&dt, "(main) hash error")
   ;  }

      if (INTERFACE(par, MODE_STOP_WORDS))
      {  char buf[16]
      ;  for (i=0;i<dt.wordlist.n_words;i++)
         fprintf
         (  stdout
         ,  "%s\t%d\t%.6f\n"
         ,  get_sylmer(par.K_g,  dt.wordlist.words[i], buf)
         ,  (int) dt.wordlist.words[i]
         ,  get_word_entropy(dt.wordlist.words[i], par.K_g)
         )
      ;  tell
         (  VERBOSE_ONCE
         ,  "printed %d words (K=%d)"
         ,  (int) dt.wordlist.n_words
         ,  (int) par.K_g
         )
      ;  return myfinish(&dt, 0)
   ;  }



   /* ***************************
  *   more word stuff above.
 */

      if (par.gap_min <= par.gap_max)
      par.modes |= MODE_GAPPY

   ;  if (par.M_g && !(INTERFACE(par, MODE_R2USERON)))
      par.modes |= MODE_R2ON, par.modes ^= MODE_R2ON

   ;  if (INTERFACE(par, MODE_GAPPY))
      {  if (par.K_g & 1)
         return mainfail(&dt, "gapped mode requires -k parameter to be even")
      ;  if
         (  !(par.dyad_modes & DYAD_UNIFORM)
         && par.M_g
         && par.M_g * 2 > par.K_g
         )
            miaow("in gap mode without --dyad-uniform -m parameter must not exceed K/2, (fixed it for you)")
         ,  par.M_g = par.K_g / 2
      ;  else if (par.dyad_modes & DYAD_UNIFORM)
         par.M_g = par.K_g / 2
   ;  }

      if (!(par.fpresult = myfopen(par.fnresult, "w")))
      return mainfail(&dt, "unable to express my sincerest P-values")

   ;  if (MAIN_MODE(par, MAIN_COUNT | MAIN_WSTAT) && INTERFACE(par, MODE_R2ON))
      miaow("--r2-off not set in word mode")

   ;  if (INTERFACE(par, MODE_LOGFOLD) && !MAIN_MODE(par, MAIN_TABLE))
         miaow("--logfold implies --table mode (now set)")
      ,  par.main_modes |= MAIN_TABLE

   ;  if (!INTERFACE(par, MODE_OVER | MODE_UNDER | MODE_NONE | MODE_LOGFOLD))
      par.modes |= MODE_OVER | MODE_UNDER

   ;  {  const char* wrong = NULL

      ;  if (par.splitsize && par.binsize)
         wrong = "conflicting bin options (-grow/step & -grow/step-times)"

      ;  if (par.fnexp && INTERFACE(par, MODE_FREQUENCIES))
         wrong = "cannot combine -dump-expected and -expected"

      ;  if (par.K_g > 15)
         wrong = "-k option not in range [0,15]"

      ;  if (par.M_g && (par.M_g < 1 || par.M_g > par.K_g))
         wrong = "-m value not in range [1, K]"

      ;  if (INTERFACE(par, MODE_FREQUENCIES) && par.M_g)
         wrong = "-expected and -m do not jive"

      ;  if (par.M_g && (par.per_sequence || par.cap))
         wrong = "(-threshold/-cap) -m do not shimmy"

      ;  if (INTERFACE(par, MODE_BYSTRETCH) && par.per_sequence)
         wrong = "--ck and -threshold do not walz"

      ;  if (MAIN_MODE(par, MAIN_EXTREME) && MAIN_MODE(par, MAIN_TABLE))
         par.main_modes ^= MAIN_TABLE

      ;  if (wrong)
         return mainfail(&dt, wrong)
   ;  }

      if (MAIN_MODE(par, MAIN_WSTAT))
      {  if (!par.M_g)
         par.M_g = par.K_g == 1 ? 1 : par.K_g - 1
      ;  if (0 && !par.M_g)      /* recentlychanged */
         par.M_g = 1
   ;  }

/* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */

      if (par.fpquery)
      {  status = read_query(par.fpquery, dt.thequery, INTERFACE(par, MODE_REVERSE))
      ;  myfzclose(&par.fpquery)

      ;  if (status != SYL_DONE)
         return mainfail(&dt, "no entries found in query file")

         /* we need to sort query so that we bsearch it in read_fasta,
          * as well as for making sure that query IDs are unique.
          * 'encodes size of query' attribute is moved from query_n_file to
          * query_n_uniq below. that's somewhat confusing, but we need both.
          * ownership of query is not clearly localized (namely here, with
          * the sorting and whatnot, and in query.c). That's Bad.
         */
      ;  query_uniq(dt.thequery)
      ;  if (dt.thequery->qel_n_uniq < dt.thequery->qel_n_file)
         miaow
         (  "%lu duplicate entries in query file"
         ,  (unsigned long) (dt.thequery->qel_n_file - dt.thequery->qel_n_uniq)
         )
               /* _ we could plough on, but bsearch in read_fasta will crash.
                * Let's not pursue this corner case. Fixme: stale comment or code.
               */
      ;  if (dt.thequery->qel_n_uniq == 0)
         {
         }
   ;  }

      status
      =  read_fasta
         (  dt.thecollection
         ,  par.fpfasta
         ,  dt.thequery
         ,  &dt.rank_max
         ,  par.modes
         ,  par.length_co
         ,  par.repeat_mask
         ,  INTERFACE(par, MODE_GAPPY) ? 1 + (par.K_g >> 1) : 0
         ,  INTERFACE(par, MODE_GAPPY) ? par.K_g >> 1 : par.K_g
         )
   ;  myfzclose(&par.fpfasta)

   ;  if (status)
      return mainfail(&dt, "failed reading fasta file")

               /* Can we continue as a stress test?
                * divide by zero lurks in too many places.
               */
   ;  if (!dt.thecollection->list_n)
      {  exit_g = 4
      ;  return mainfail
         (  &dt
         ,  "nothing left in fasta file (%s mode)"
         ,  INTERFACE(par, MODE_UNIVERSE) ? "universe" : "subset"
         )
   ;  }

      if (INTERFACE(par, MODE_STOP_READ))
      return myfinish(&dt, 0)

                  /* write missing and/or found entries if requested to
                   * hrm. we qsort below, duplicated further down
                   * fixme Move this to appropriate module
                  */
   ;  if ((par.fnmissing || par.fnfound || par.fnwquery) && dt.thequery->qel_list)
      {  FILE* iffpmissing = par.fnmissing ? myfopen(par.fnmissing, "w") : NULL
      ;  FILE* iffpfound = par.fnfound ? myfopen(par.fnfound, "w") : NULL
      ;  FILE* iffpqueryout = par.fnwquery ? myfopen(par.fnwquery, "w") : NULL
      ;  qsort
         (  dt.thequery->qel_list
         ,  dt.thequery->qel_n_uniq
         ,  sizeof dt.thequery->qel_list[0]
         ,  qel_cmp_by_rank
         )
      ;  for (i=0;i<dt.thequery->qel_n_uniq;i++)
         {  if (dt.thequery->qel_list[i].found && iffpfound)
            fprintf(iffpfound, "%s\n", dt.thequery->qel_list[i].sid)
         ;  else if (!dt.thequery->qel_list[i].found && iffpmissing)
            fprintf(iffpmissing, "%s\n", dt.thequery->qel_list[i].sid)
         ;  if (iffpqueryout)
            fprintf(iffpqueryout, "%s\n", dt.thequery->qel_list[i].sid)
      ;  }
         if (iffpmissing)    myfclose(&iffpmissing)
      ;  if (iffpfound)      myfclose(&iffpfound)
      ;  if (iffpqueryout)   myfclose(&iffpqueryout)
   ;  }

      if (dt.thequery->qel_list && !dt.thequery->qel_n_found)
      {  exit_g = 2
      ;  return mainfail(&dt, "no query IDs found in fasta file - I give up")
   ;  }

      if (!par.K_g)
      return 0

         /* This doubles the collection of sequences, with input sequences and
          * their reverse complement interweaved. Any code that wants to
          * distinguish between the two types can do so by looking at the
          * parity of the index. This is a simple-minded approach that retains
          * all information and allows maximum flexibility.
          * Needs a contract wrt identifiers and sorting though.
         */
   ;  if (INTERFACE(par, MODE_2STRAND))
      {  if (collection_revcompl(dt.thecollection))
         return mainfail(&dt, NULL)
      ;  dt.collection_n_found = 2 * dt.thequery->qel_n_found
   ;  }
      else if (INTERFACE(par, MODE_2STRAND2))
      {  if (collection_revcompl2(dt.thecollection, par.gap_max))
         return mainfail(&dt, NULL)
      ;  dt.collection_n_found = dt.thequery->qel_n_found
   ;  }
      else
      dt.collection_n_found = dt.thequery->qel_n_found


                   /* now put query elements in front in the right order
                    * query reversal is already done above
                    * fixme: put in module.
                    * TODO: better encapsulation/ownership for the main
                    * constants, such as dt.collection_n_found and
                    * dt.thecollection->list_n
                   */
   ;  if (dt.thequery->qel_list)
      {  qsort
         (  dt.thecollection->list
         ,  dt.thecollection->list_n
         ,  sizeof dt.thecollection->list[0]
         ,  unit_cmp_by_rank
         )
                  /* if asked get rid of ulterior elements */
      ;  if (par.shift_left+par.shift_right)
         {  if (par.shift_left + par.shift_right >= dt.collection_n_found)
            miaow("shift exceeds number query elements found!")
         ;  else
            {  for (i=0;i<par.shift_left;i++)
               dt.thecollection->list[i].rank = dt.rank_max+i
            ;  for (i=0;i<par.shift_right;i++)
               dt.thecollection->list[dt.thecollection->list_n-par.shift_right+i].rank = dt.rank_max+par.shift_left+i
            ;  qsort                      /* fixme: inefficient but oh well */
               (  dt.thecollection->list
               ,  dt.thecollection->list_n
               ,  sizeof dt.thecollection->list[0]
               ,  unit_cmp_by_rank
               )
            ;  dt.collection_n_found    -= (par.shift_left+par.shift_right)
            ;  dt.thecollection->list_n -= (par.shift_left+par.shift_right)
                                          /* fixme: dt.rank_max now meaningless */
         ;  }
            tell
            (  VERBOSE_ONCE
            ,  "shifted sequences/queries now %lu/%lu"
            ,  (unsigned long) dt.thecollection->list_n
            ,  (unsigned long) dt.collection_n_found
            )
      ;  }
      }
      else  /* no query */
      {  if (INTERFACE(par, MODE_SORTBYLENGTH))
         qsort
         (  dt.thecollection->list
         ,  dt.thecollection->list_n
         ,  sizeof dt.thecollection->list[0]
         ,  INTERFACE(par, MODE_REVERSE) ? unit_cmp_by_revlength : unit_cmp_by_length
         )
      ;  else if (INTERFACE(par, MODE_REVERSE))
         qsort
         (  dt.thecollection->list
         ,  dt.thecollection->list_n
         ,  sizeof dt.thecollection->list[0]
         ,  unit_cmp_by_revrank
         )
   ;  }


/* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */

            /* fixme. BTB, block too big.
             * mitigating: several modes that need shared functionality,
             * allowing many of the same parameters and options (but not all).
             * It's cumbersome and big and should be improved, but it is not a
             * complete disaster.  There are some dependencies, e.g.
             * dt.wordlist.n_words is used in computing dt.table dimensions.
            */

      {  dim n_universe = dt.thecollection->list_n
      ;  dim n_bin_g = 0
      ;  dim subset_size      =  0
      ;  sylstatus status = SYL_SUCCESS

      ;  mkv* markov_b, *markov_l   /* Background/gloBal, LocaL */

      ;  acc threshold              /* manage per sequence thresholded counts */

      ;  acc* acc_threshold = par.per_sequence || par.cap ? &threshold : NULL

      ;  dim deus_ex_machina = 0    /* times we have to readjust expected occurrence */
      ;  dim sudo_ex_machina = 0    /* times we have to readjust because of par.unit_size_fake */

      ;  dim bg_count = 0
      ;  dim x

      ;  dim n_analyse = dt.thequery->qel_list ? dt.collection_n_found : n_universe

      ;  if
         (  MAIN_MODE(par, MAIN_EXTREME)
         && !(dt.extrema = mymalloc(dt.wordlist.n_words * sizeof dt.extrema[0] * 2))
         )
         return mainfail(&dt, "could not alloc extrema")

      ;  if (INTERFACE(par, MODE_UNIVERSE))
         n_universe = n_analyse
                                    /* n_analyse is the exact number of sequences that
                                     * we analyze, bin-wise.
                                     *
                                     * n_universe is the size of the background.
                                     * Unless -subset is used, it will be the same as n_analyse.
                                     *
                                     * dt.collection_n_found is the exact set-size of those query identifiers
                                     * that were found in the background file.
                                    */
      ;  if (par.M_g)
            markov_b = &dt.mkv1
         ,  markov_l = &dt.mkv2
      ;  else
            markov_b = NULL
         ,  markov_l = NULL

      ;  if (!n_analyse)
         return mainfail(&dt, "no units, no party")
      ;  dt.n_analyse = n_analyse
      ;  dt.n_universe = n_universe

                     /* hierverder: require par.binsize for _mww */
      ;  if (par.splitsize || par.binsize)
         {  n_bin_g
            =  interval_prepare
               (  dt.thecollection
               ,  n_analyse
               ,  par.splitsize
               ,  par.binsize
               ,  INTERFACE2(par, Mode_BINBYLENGTH)
               )
                  /* fixme: this condition hurts the eye. collate some stuff somewhere */
         ;  if
            (  n_bin_g == 1
            && dt.thequery->qel_list
            && n_analyse == n_universe
            && (  !par.n_universe_fake
               && !par.utimes_fake
               && !MAIN_MODE(par, MAIN_COUNT)
               && !par.format
               )
            )
            miaow("___[warning] single bin on entire universe renders hypergeometric useless")
      ;  }

         else if (INTERFACE(par, MODE_FSTAT))
         {  acc ctwindow
         ;  if (acc_init(&ctwindow, par.K_g, dt.hash_template_g))
            return mainfail(&dt, NULL)
         ;  for (i=0;i<n_universe;i++)
            {  unit* ut = dt.thecollection->list+i
            ;  acc_addto(par.K_g, ut, &ctwindow, 0)
         ;  }
            tell
            (  VERBOSE_ONCE
            ,  "bases %lu, masked %.1f, masked k-mers %.1f"
            ,  (unsigned long) ctwindow.n_bases
            ,  (double) (ctwindow.n_mask * 100.0 / (ctwindow.n_bases + 1))
            ,  (double) (ctwindow.n_unsyl * 100.0 / (ctwindow.n_syl + ctwindow.n_unsyl))
            )
         ;  return myfinish(&dt, 0)
      ;  }
         else if (!MAIN_MODE(par, MAIN_COUNT | MAIN_FILTER | MAIN_WSTAT | MAIN_MWW) && !par.format)
         return mainfail(&dt, "I do not know how to make bins")

      ;  dt.n_bin = n_bin_g

                        /* fixme; clause below is unwieldy. organize. */
      ;  if
         (  MAIN_MODE(par, MAIN_TABLE)
            && !par.format
            && !MAIN_MODE(par, MAIN_COUNT | MAIN_FILTER)
         )
         {  if (!(dt.table = mycalloc(n_bin_g * dt.wordlist.n_words, sizeof dt.table[0])))
            return mainfail(&dt, "could not alloc table")
         ;  if (!(dt.bins_info = mycalloc(n_bin_g, sizeof dt.bins_info[0])))
            return mainfail(&dt, "could  not alloc header!")
      ;  }

         if (acc_init(&dt.bg, par.K_g, dt.hash_template_g))
         return mainfail(&dt, NULL)
      ;  if (acc_init(&dt.window, par.K_g, dt.hash_template_g))
         return mainfail(&dt, NULL)
      ;  if (INTERFACE(par, MODE_GAPPY) && acc_init(&dt.bg_gap_cache, par.K_g, dt.hash_template_g))
         return mainfail(&dt, NULL)
      ;  if (INTERFACE(par, MODE_GAPPY) && acc_init(&dt.window_gap_cache, par.K_g, dt.hash_template_g))
         return mainfail(&dt, NULL)

      ;  if (acc_threshold && acc_init(acc_threshold, par.K_g, dt.hash_template_g))
         return mainfail(&dt, NULL)

      ;  if (par.M_g)
         {  if (acc_init(&dt.bg_xy ,par.M_g-1, NULL)) return mainfail(&dt, NULL)
         ;  if (acc_init(&dt.bg_xyz,  par.M_g, NULL)) return mainfail(&dt, NULL)
         ;  if (acc_init(&dt.window_xy ,par.M_g-1, NULL)) return mainfail(&dt, NULL)
         ;  if (acc_init(&dt.window_xyz,  par.M_g, NULL)) return mainfail(&dt, NULL)
         ;  if (mkv_init(markov_l, par.K_g, dt.hash_template_g)) return mainfail(&dt, NULL)
         ;  if (mkv_init(markov_b, par.K_g, dt.hash_template_g)) return mainfail(&dt, NULL)
      ;  }

                  /* Sequence mode.
                   * should be split off into different program.
                  */
         if (par.format)
         {  if (strstr(par.format, "%C"))      /* fixme this errs mistakenly on "%%C" */
            miaow("___[sequence] possibly useless use of %%C directive (taking constant 0)")
         ;  for (i=0;i<n_universe;i++)
            {  if (par.lastbin && i >= par.lastbin)
               break
            ;  unit_filter(par.fpresult, par.format, dt.thecollection->list+i, NULL, 0)
         ;  }
            return myfinish(&dt, 0)
      ;  }

                  /* Filter mode illustrates the basic loop.
                   * P-value mode, further below, is monstrous in comparison.
                   * NOTE: does/should not use bg_count.
                  */
         else if (MAIN_MODE(par, MAIN_COUNT | MAIN_FILTER))
         return dispatch_count(&par, &dt)

                  /* Experimental
                  */
      ;  else if (MAIN_MODE(par, MAIN_MWW))
         {  if (dispatch_mww(&par, &dt, INTERFACE(par, MODE_WSTAT2) ? 1 : 0))
            return mainfail(&dt, "mww mode failed")
         ;  else
            return myfinish(&dt, 0)
      ;  }

                  /* This block accumulates count for the background.
                   * It also detects the subset size subset_size.
                   *
                   * fixme: subset_size logic. (document) what happens when
                   * subset_size == 0?
                  */
       ; for (x=0;x<n_universe;x++)
         {  unit* ut = dt.thecollection->list+x

         ;  if (par.per_sequence || par.cap)
            {  if (acc_addto(par.K_g, ut, acc_threshold, INTERFACE(par, MODE_R2ON)))
               return mainfail(&dt, "addto threshold background")
            ;  if (acc_per_sequence(par.K_g, &dt.bg, acc_threshold, ut, par.per_sequence, par.cap))
               return mainfail(&dt, "perseq background")
         ;  }
            else
            {  status
               =  INTERFACE(par, MODE_GAPPY)
               ?  acc_addto_gapped(par.K_g, ut, &dt.bg, par.gap_min, par.gap_max, par.gap_gap, &dt.bg_gap_cache)
               :  acc_addto(par.K_g, ut, &dt.bg, INTERFACE(par, MODE_R2ON))
            ;  if (status)
               return mainfail(&dt, "addto background")
         ;  }

            if (par.M_g)
            {  if (acc_addto(par.M_g,   ut, &dt.bg_xyz, INTERFACE(par, MODE_R2ON)))
               return mainfail(&dt, "addto xyz background")
            ;  if (acc_addto(par.M_g-1, ut, &dt.bg_xy , INTERFACE(par, MODE_R2ON)))
               return mainfail(&dt, "addto xy background")
         ;  }

            if (x+1 == dt.collection_n_found)
            {  if (INTERFACE(par, MODE_BYSTRETCH))
               subset_size = dt.bg.n_spaced
            ;  else if (par.unit_size_fake)
               subset_size = dt.bg.n_unit * par.unit_size_fake
            ;  else if (par.per_sequence)
               subset_size = dt.bg.n_unit
            ;  else
               subset_size = dt.bg.n_syl
         ;  }
         }

         if (par.n_universe_fake)
            bg_count = par.n_universe_fake
      ;  else if (par.utimes_fake * subset_size)
            bg_count = par.utimes_fake * subset_size

/* the four options below
 * also apply to the window size
*/
      ;  else if (par.unit_size_fake)
            bg_count = n_universe * par.unit_size_fake
      ;  else if (par.per_sequence)
            bg_count = dt.bg.n_unit
      ;  else if (INTERFACE(par, MODE_BYSTRETCH))
            bg_count = dt.bg.n_spaced
      ;  else
            bg_count = dt.bg.n_syl

      ;  if (par.n_universe_fake || par.utimes_fake)
         tell
         (  VERBOSE_ONCE
         ,  "background count %lu imposed count %lu"
         ,  (unsigned long) dt.bg.n_syl
         ,  (unsigned long) bg_count
         )

      ;  if (!bg_count || !n_universe)
         return mainfail(&dt, "no background count or no universe")

      ;  dt.bg_use_count = bg_count

                     /* output the universe counts as frequencies
                      * fixme funcify
                      * make sure works with -subset/-universe/-do.
                      * use par.fpresult
                      * This is a nice small example of how things work.
                     */
      ;  if (MAIN_MODE(par, MAIN_WSTAT))
         return dispatch_wstat(&par, &dt)

      ;  if (verbose_G)
         {  dim sylsum = acc_get_sum(&dt.bg)

         ;  tell
            (  VERBOSE_ONCE
            ,  "bases %lu, masked %lu, counted %lu, rejected %lu"
               ", spaced %lu, r2-skipped %lu, universe %lu"
            ,  (unsigned long) dt.bg.n_bases
            ,  (unsigned long) dt.bg.n_mask
            ,  (unsigned long) sylsum
            ,  (unsigned long) dt.bg.n_unsyl
            ,  (unsigned long) dt.bg.n_spaced
            ,  (unsigned long) dt.bg.n_r2
            ,  (unsigned long) dt.bg.n_syl
            )
      ;  }

         if (par.M_g)
         {  status 
            =  INTERFACE(par, MODE_GAPPY)
            ?   
               have_expectations_gapped
               (  par.fnexp
               ,  0                 /* signifies background */
               ,  par.K_g
               ,  par.M_g
               ,  dt.bg.n_syl       /* fixme quantify edge effects (it's not consecutive n_syl) */
               ,  par.gap_max - par.gap_min + 1
               ,  dt.bg_xy.acount->cells
               ,  dt.bg_xyz.acount->cells
               ,  markov_b
               ,  par.dyad_modes & DYAD_UNIFORM ? 0 : 1
               )
            :
               have_expectations
               (  par.fnexp
               ,  0                 /* signifies background */
               ,  par.K_g
               ,  par.M_g
               ,  dt.bg_xy.acount->cells
               ,  dt.bg_xyz.acount->cells
               ,  markov_b
               )
         ;  if (status)
            return mainfail(&dt, "failed to compute background expectations")
      ;  }

         if (INTERFACE(par, MODE_STOP_BG))
         return myfinish(&dt, 0)


/* If trial == 1 we have a regular sylamer run.
 * This is the main monster loop that could do
 * with some folding. monsterloop!
*/
      ;  {  for (trial_i = 0; trial_i < par.n_trial; trial_i++)
            {  dim window_count = 0
            ;  dim window_count_prev = 0
            ;  char sylbuf[16]
            ;  spacer sp = { 0 }

            ;  if (!trial_i && !MAIN_MODE(par, MAIN_TABLE))
               fprintf
               (  par.fpresult
               ,  "binid\tword\tsgnlogh\tsgnlogb\tn1\tn1n2\tk\tndrawn"
                  "\tfreqbg\tfreqwdw\tratio\tlogro\tn1exp%s\n"
               ,  (INTERFACE(par, MODE_PRINTID)) ? "\tseqid" : ""
               )

            ;  if (trial_i)
               {  acc_reset(&dt.window)
               ;  if (par.M_g)
                  {  acc_reset(&dt.window_xy)
                  ;  acc_reset(&dt.window_xyz)
               ;  }
               }
                     /* shuffle unless:
                      * -  there is more than one trial,
                      * -  we have a query (universe/subset) file
                      * -  we are at the first trial
                     */
               if (INTERFACE(par, MODE_SHUFFLE) && !(dt.thequery->qel_list && par.n_trial > 1 && !trial_i))
               collection_shuffle(dt.thecollection, n_universe)

            ;  if (dt.extrema)
               for (x=0;x<dt.wordlist.n_words;x++)
                  dt.extrema[2*x] = 1.0
               ,  dt.extrema[2*x+1] = 1.0

            ;  for (x=0;x<n_bin_g;x++)       /* the bin-loop */
               {  acc_iter iterator
               ;  dim i

               ;  if (par.lastbin && x >= par.lastbin)
                  break

                                             /* compute i_lo and i_hi */
               ;  interval_step(dt.thecollection, dt.n_analyse, &sp)
;if(0)fprintf(stderr, "invoking interval_step\n")

               ;  if (!trial_i && dt.bins_info)
                  {  dt.bins_info[x].offset = sp.i_hi
                  ;  dt.bins_info[x].size   = sp.i_hi - sp.i_lo
                  ;  dt.bins_info[x].utr_length = 0.0
                  ;  dt.bins_info[x].n_masked   = 0.0
                  ;  dt.bins_info[x].gccontent  = 0.0
               ;  }

                                             /* reset all the counts for regular bins */
                  if (x && INTERFACE2(par, Mode_STEPPING))
                  {  acc_reset(&dt.window)
                  ;  if (par.M_g)
                     {  acc_reset(&dt.window_xy)
                     ;  acc_reset(&dt.window_xyz)
                  ;  }
                  }
                                             /* update the counts */
                  for (i=sp.i_lo;i<sp.i_hi;i++)
                  {  unit* ut = dt.thecollection->list+i

                  ;  if (!trial_i && dt.bins_info)
                     {  dt.bins_info[x].utr_length += ut->n_bases
                     ;  dt.bins_info[x].n_masked   += ut->n_masked
                     ;  dt.bins_info[x].gccontent  += ut->gccontent
                  ;  }

                     if (par.per_sequence || par.cap)
                     {  if (acc_addto(par.K_g, ut, acc_threshold, INTERFACE(par, MODE_R2ON)))
                        return mainfail(&dt, "addto threshold window")
                     ;  if (acc_per_sequence(par.K_g, &dt.window, acc_threshold, ut, par.per_sequence, par.cap))
                        return mainfail(&dt, "perseq window")
                  ;  }
                     else
                     {  status
                        =  INTERFACE(par, MODE_GAPPY)
                        ?  acc_addto_gapped(par.K_g, ut, &dt.window, par.gap_min, par.gap_max, par.gap_gap, &dt.window_gap_cache)
                        :  acc_addto(par.K_g, ut, &dt.window, INTERFACE(par, MODE_R2ON))
                     ;  if (status)
                        return mainfail(&dt, "addto window")
                  ;  }
                              /* fixme: make sure gapped mode and Markov is predictable and reliable */
                     if (par.M_g)
                     {  if (acc_addto(par.M_g,    ut, &dt.window_xyz, INTERFACE(par, MODE_R2ON)))
                        return mainfail(&dt, "addto xyz window")
                     ;  if (acc_addto(par.M_g-1,  ut, &dt.window_xy , INTERFACE(par, MODE_R2ON)))
                        return mainfail(&dt, "addto xy window")
                  ;  }
                  }

                  if (par.unit_size_fake)
                     window_count = dt.window.n_unit * par.unit_size_fake
               ;  else if (par.per_sequence)
                     window_count = dt.window.n_unit
               ;  else if (INTERFACE(par, MODE_BYSTRETCH))
                     window_count = dt.window.n_spaced
               ;  else
                     window_count = dt.window.n_syl

               ;  if (verbose_G)
                  {  if (verbose_G == VERBOSE_ONCE)
                     fputc('.', stderr)
                  ;  else if (INTERFACE2(par, Mode_STEPPING))
                     tell
                     (  VERBOSE_BIN
                     ,  "bin %lu [%lu..%lu) size %lu/%lu (%.2f)"
                     ,  (unsigned long) (x+1)
                     ,  (unsigned long) sp.i_lo
                     ,  (unsigned long) sp.i_hi
                     ,  (unsigned long) (sp.i_hi - sp.i_lo)
                     ,  (unsigned long) window_count
                     ,  (double) (100.0f * window_count) / (1.0 * bg_count)
                     )
                  ;  else
                     tell
                     (  VERBOSE_BIN
                     ,  "bin %lu [0..%lu) adding %lu/%lu, total %lu (%.2f)"
                     ,  (unsigned long) (x+1)
                     ,  (unsigned long) sp.i_hi
                     ,  (unsigned long) (sp.i_hi - sp.i_hi_prev)
                     ,  (unsigned long) (window_count - window_count_prev)
                     ,  (unsigned long) window_count
                     ,  (double) (100.0f * window_count) / (1.0 * bg_count)
                     )
               ;  }

                  window_count_prev = window_count
               ;  sp.i_hi_prev = sp.i_hi

               ;  if (par.M_g)
                  {  status
                     =  INTERFACE(par, MODE_GAPPY)
                     ?
                        have_expectations_gapped
                        (  n_bin_g == 1 ? NULL : par.fnexp
                        ,  x+1
                        ,  par.K_g
                        ,  par.M_g
                        ,  dt.window.n_syl
                        ,  par.gap_max - par.gap_min + 1
                        ,  dt.window_xy.acount->cells
                        ,  dt.window_xyz.acount->cells
                        ,  markov_l
                        ,  par.dyad_modes & DYAD_UNIFORM ? 0 : 1
                        )
                     :
                        have_expectations
                        (  n_bin_g == 1 ? NULL : par.fnexp
                        ,  x+1
                        ,  par.K_g
                        ,  par.M_g
                        ,  dt.window_xy.acount->cells
                        ,  dt.window_xyz.acount->cells
                        ,  markov_l
                        )
                  ;  if (status)
                     return mainfail(&dt, "failed to compute window %d expectations", (int) x)
               ;  }

                  acc_iter_reset(&iterator)    /* we're at bin-loop level */

                     /* the iterator abstracts over that we either
                      * have a regular array or a hash array. It may throw an error in &status
                      * if its internal computations do not add up - that would be a bug, that would.
                     */
               ;  while (acc_bi_iterate(&iterator, &dt.bg, &dt.window, markov_b, markov_l, dt.eoracle_g, &status))
                  {  dim bg_count_y    =  iterator.count1
                  ;  dim window_count_y=  iterator.count2
                  ;  float eglobal_y   =  iterator.val1
                  ;  float elocal_y    =  iterator.val2
                  ;  float eoracle_y   =  iterator.exp

                  ;  dim y =  iterator.hid            /* should never be used as offset now */
                  ;  double hyperQ = 1.0, hyperP = 1.0, binomQ = 1.0, binomP = 1.0

                                       /* some ugliness for markov pretended background count */
                  ;  double anchor_fac
                     =     par.M_g && !INTERFACE(par, MODE_NOANCHOR) && eglobal_y
                        ?     bg_count_y * 1.0
                           / (1.0 * eglobal_y * bg_count)
                        :  1.0

                  ;  double P =     markov_l
                                 ?  elocal_y * anchor_fac
                                 :     dt.eoracle_g
                                    ?  eoracle_y
                                    :  -1.0

                  ;  dim bg_count_yy
                     =     par.M_g
                        ?  elocal_y * bg_count * anchor_fac
                        :     dt.eoracle_g
                           ?  eoracle_y * bg_count
                           :  bg_count_y

                  ;  dim bg_count_yy_exp = bg_count_yy

                  ;  double fratio = bg_count_yy && window_count
                     ?  (1.0 * window_count_y / window_count) / (1.0 * bg_count_yy / bg_count)
                     :  1.0

                  ;  double flogratio = fratio ? log(fratio) / log(2.0) : -10.0

                  ;  if (P > 1.0)
                        miaow("big P %g for word %s", P, get_sylmer(par.K_g, y, sylbuf))
                     ,  P = 1.0

                                       /* This may happen, anchored markov is not
                                        * guarantueed to exceed window_count_y
                                       */
                  ;  if (bg_count_yy < window_count_y)
                     {  if (par.iffplog)
                        fprintf
                        (  par.iffplog
                        ,  "deus %s %lu %lu %lu\n"
                        ,  get_sylmer(par.K_g, y, sylbuf)
                        ,  (unsigned long) (x+1)
                        ,  (unsigned long) bg_count_yy
                        ,  (unsigned long) window_count_y
                        )
                     ;  bg_count_yy = window_count_y
                     ;  deus_ex_machina++
                  ;  }
                                       /* This may happen e.g. with par.unit_size_fake */
                     if (window_count < window_count_y || bg_count < bg_count_yy)
                     {  if (par.iffplog)
                        fprintf
                        (  par.iffplog
                        ,  "sudo %s %lu w=%lu wy=%lu b=%lu by=%lu el=%g eg=%g\n"
                        ,  get_sylmer(par.K_g, y, sylbuf)
                        ,  (unsigned long) (x+1)
                        ,  (unsigned long) window_count
                        ,  (unsigned long) window_count_y
                        ,  (unsigned long) bg_count
                        ,  (unsigned long) bg_count_yy
                        ,  elocal_y
                        ,  eglobal_y
                        )
                     ;  sudo_ex_machina++

                     ;  if (window_count < window_count_y)
                        window_count_y = window_count
                     ;  if (bg_count < bg_count_yy)
                        bg_count_yy = bg_count
                  ;  }

#ifdef ARRR
#  define LOWER_TAIL 1
#  define UPPER_TAIL 0
#  define NO_LOG     0
#  define   get_hyper_q(white_sample, white_urn, black_urn, sample)  \
               phyper(white_sample, white_urn, black_urn, sample, UPPER_TAIL, NO_LOG)
#  define   get_hyper_p(white_sample, white_urn, black_urn, sample)  \
               phyper(white_sample, white_urn, black_urn, sample, LOWER_TAIL, NO_LOG)
#  define   get_binom_q(white_sample, P, sample)   pbinom(white_sample, sample, P, UPPER_TAIL, N0_LOG)
#  define   get_binom_q(white_sample, P, sample)   pbinom(white_sample, sample, P, LOWER_TAIL, N0_LOG)
#else
#  define   get_hyper_q(white_sample, white_urn, black_urn, sample)  \
               gsl_cdf_hypergeometric_Q(white_sample, white_urn, black_urn, sample)
#  define   get_hyper_p(white_sample, white_urn, black_urn, sample)  \
               gsl_cdf_hypergeometric_P(white_sample, white_urn, black_urn, sample)
#  define   get_binom_q(white_sample, P, sample)   gsl_cdf_binomial_Q(white_sample, P, sample)
#  define   get_binom_p(white_sample, P, sample)   gsl_cdf_binomial_P(white_sample, P, sample)
#endif

;if (0 && x+1 == 10)
fprintf(stdout, "   bg exp (%d) win ct (%d) win size (%d) fac (%.2f) base fac (%.4f)\n",
(int) bg_count_yy, (int) window_count_y, (int) window_count, anchor_fac, (double) window_count * 1.0 / bg_count)
;
                     if (INTERFACE(par, MODE_OVER) && window_count_y)
                        hyperQ
                        =  get_hyper_q
                           (  window_count_y-1           /* > k-1 in window  */
                           ,  bg_count_yy                /* n1 in population */
                           ,  bg_count - bg_count_yy     /* n2 in population */
                           ,  window_count               /* window size      */
                           )

                  ;  if (INTERFACE(par, MODE_OVER) && !MAIN_MODE(par, MAIN_TABLE) && INTERFACE(par, MODE_BINOMIAL) && window_count_y)
                        binomQ
                        =  get_binom_q
                           (  window_count_y-1
                           ,  P >= 0 ? P : bg_count_yy * 1.0 / bg_count
                           ,  window_count
                           )

                  ;  if (INTERFACE(par, MODE_UNDER) && window_count_y < window_count)
                        hyperP
                        =  get_hyper_p
                           (  window_count_y             /* <= k in window   */
                           ,  bg_count_yy                /* n1 in population */
                           ,  bg_count - bg_count_yy     /* n2 in population */
                           ,  window_count               /* window size      */
                           )

;if(0)fprintf(stderr, "P= %.5g, elocal_y = %.5g\n", (double) P,  (double) elocal_y)
                  ;  if (INTERFACE(par, MODE_UNDER) && !MAIN_MODE(par, MAIN_TABLE) && INTERFACE(par, MODE_BINOMIAL) && window_count_y < window_count)
                        binomP
                        =  get_binom_p
                           (  window_count_y
                           ,  P >= 0 ? P : bg_count_yy * 1.0 / bg_count
                           ,  window_count
                           )
;if(0)fprintf(stdout, "P=%.5g bgcountf=%.3g bgcount=%u, hyper=%.5g,%.5g binom=%.5g,%.5g\n", (double) P, (double) (elocal_y * bg_count * anchor_fac), (unsigned) bg_count_yy, (double) hyperQ, (double) hyperP, (double) binomQ, (double) binomP)
;if(0)fprintf(stderr, "bgcount = %u, hyper = %.5g, %.5g binom = %.5g, %.5g\n", (unsigned) bg_count_yy, (double) hyperQ, hyperP, binomQ, binomP)

                  ;  if (dt.extrema)
                     {  dim offset = 2 * (iterator.offset-1)
                     ;  if (hyperP < dt.extrema[offset])
                        dt.extrema[offset] = hyperP
                     ;  if (hyperQ < dt.extrema[offset+1])
                        dt.extrema[offset+1] = hyperQ
                  ;  }
                     else if (MAIN_MODE(par, MAIN_TABLE))
                     {  double h = hyperP < hyperQ ? hyperP : hyperQ
                     ;  if (INTERFACE(par, MODE_LOGFOLD))
                        h = flogratio
                     ;  else
                        {  double h2 = h > 0 ? - (log(h) / log(10.0)) : out_of_bounds_g
                        ;  if (h2 > out_of_bounds_g)
                           h2 = out_of_bounds_g
                        ;  if (  (!h && 2 * bg_count_yy > bg_count)
                              || hyperP < hyperQ
                              )
                           h2 = -h2
                        ;  h = h2
                     ;  }
                        dt.table[(iterator.offset-1) * n_bin_g + x] = h
                  ;  }
                     else if (!par.cutoff || hyperP < par.cutoff || hyperQ < par.cutoff)
                     {  double h = hyperQ < hyperP ? hyperQ : hyperP
                     ;  double b = binomQ < binomP ? binomQ : binomP      /* instead, use same test for consistency? */
                     ;  double normalise_factor = par.per_sequence ? 1.0 : LSIZE(par.K_g)

                     ;  h = log(h) / log(10.0)
                     ;  b = log(b) / log(10.0)

                     ;  if (hyperQ < hyperP)      /* it's overrepresented */
                        h *= -1.0

                     ;  if (binomQ < binomP)      /* it's overrepresented */
                        b *= -1.0
                        
                     ;  fprintf
                        (  par.fpresult
                        ,  "%lu\t%s\t"  "%g\t%g\t" "%lu\t%lu\t%lu\t%lu\t" "%.4f\t%.4f\t%.4f\t%.4f" "\t%lu"
                        ,  (unsigned long) (sp.i_hi)
                        ,  get_sylmer(par.K_g, y, sylbuf)

                        ,  h
                        ,  b

                        ,  (unsigned long) bg_count_yy
                        ,  (unsigned long) bg_count
                        ,  (unsigned long) window_count_y
                        ,  (unsigned long) window_count

                        ,  bg_count ? normalise_factor * (bg_count_yy * 1.0 / bg_count) : 0.0
                        ,  window_count ? normalise_factor * (window_count_y * 1.0 / window_count) : 0.0
                        ,  fratio
                        ,  flogratio
                        ,  bg_count_yy_exp
                        )

                     ;  if (INTERFACE(par, MODE_PRINTID))
                        {  dim j
                        ;  fprintf(par.fpresult, "\t%s", dt.thecollection->list[sp.i_lo].sid)
                        ;  for (j=sp.i_lo+1;j<sp.i_hi;j++)
                           fprintf(par.fpresult, ",%s", dt.thecollection->list[j].sid)
                     ;  }
                        fputc('\n', par.fpresult)
                  ;  }
                  }
                  if (status)
                  return mainfail(&dt, "bi(obl)iterate")
            ;  }

               if (verbose_G == VERBOSE_ONCE)
               fputc('\n', stderr)

            ;  if (dt.extrema)
               {  dim i
               ;  for (i=0;i<dt.wordlist.n_words;i++)
                  fprintf
                  (  par.fpresult
                  ,  "%s\t%.6g\t%.6g\n"
                  ,  get_sylmer(par.K_g, dt.wordlist.words ? dt.wordlist.words[i] : i, sylbuf)
                  ,  dt.extrema[2*i]   ? -log(dt.extrema[2*i])/log(10)   : 312
                  ,  dt.extrema[2*i+1] ?  log(dt.extrema[2*i+1])/log(10) : 312
                  )
            ;  }

               if (MAIN_MODE(par, MAIN_TABLE))
               output_table(&par, &dt, trial_i)
         ;  }

            if (deus_ex_machina)
            tell
            (  VERBOSE_ONCE
            ,  "the window count exceeded the expected background times %lu times"
            ,  (unsigned long) deus_ex_machina
            )

         ;  if (sudo_ex_machina)
            tell
            (  VERBOSE_ONCE
            ,  "the word count exceeded the window size %lu times (that's very bad)"
            ,  (unsigned long) sudo_ex_machina
            )
      ;  }

         myfclose(&par.fpresult)
      ;  myfclose(&par.iffplog)
   ;  }

      return myfinish(&dt, verbose_G)
;  }


int main
(  int argc
,  char* argv[]
)
   {  return sylmain(argc, argv)
;  }



#if 0
      ;  if (MAIN_MODE(par, MAIN_WSTAT))
         {  char sylbuf[16]
         ;  dim j
         ;  acc_iter iterator

         ;  fprintf(par.fpresult, "word\twordid\tfreqbg\tfreqxp\tfreqxpa\n")

         ;  if (!par.M_g)
            return mainfail(&dt, "-m parameter required for --word-stat")

         ;  for (j=0;j<n_analyse;j++)
            {  unit* ut = dt.thecollection->list+j
            ;  if (acc_addto(par.K_g, ut, &dt.window, INTERFACE(par, MODE_R2ON)))
               return mainfail(&dt, "addto window")
            ;  if (acc_addto(par.M_g,  ut, &dt.window_xyz, INTERFACE(par, MODE_R2ON)))
               return mainfail(&dt, "addto xyz window")
            ;  if (acc_addto(par.M_g-1, ut, &dt.window_xy , INTERFACE(par, MODE_R2ON)))
               return mainfail(&dt, "addto xy window")
         ;  }

            if
            (  have_expectations
               (  NULL
               ,  1                 /* signifies not-background, ignored otherwise */
               ,  par.K_g
               ,  par.M_g
               ,  dt.window_xy.acount->cells
               ,  dt.window_xyz.acount->cells
               ,  markov_l
            )  )
            return mainfail(&dt, "failed to compute window expectations")

         ;  if
            (  have_expectations
               (  NULL
               ,  0                 /* signifies background */
               ,  par.K_g
               ,  par.M_g
               ,  dt.bg_xy.acount->cells
               ,  dt.bg_xyz.acount->cells
               ,  markov_b
            )  )
            return mainfail(&dt, "failed to compute background expectations")

         ;  acc_iter_reset(&iterator)

         ;  while (acc_bi_iterate(&iterator, &dt.bg, &dt.window, markov_b, markov_l, 0, &status))
            {  dim bg_count_y    =  iterator.count1
            ;  float eglobal_y   =  iterator.val1
            ;  float elocal_y    =  iterator.val2

            ;  double anchor_fac =  bg_count_y * 1.0 / (1.0 * eglobal_y * bg_count)
            ;  dim y =  iterator.hid            /* should never be used as offset now */

            ;  fprintf
               (  par.fpresult
               ,  "%s\t%lu\t%.4g\t%.4g\t%.4g\n"
               ,  get_sylmer(par.K_g, y, sylbuf)
               ,  y
               ,  (double) bg_count_y * 1.0 / bg_count
               ,  (double) elocal_y
               ,  (double) elocal_y * anchor_fac
               )
         ;  }
            myfclose(&par.fpresult)
         ;  return myfinish(&dt, status)
      ;  }
#endif


#if 0
         else if (MAIN_MODE(par, MAIN_COUNT | MAIN_FILTER))
         {  acc ctwindow, ctthreshold, ctcache
         ;  acc_iter iterator
         ;  acc* ctacc_threshold = par.per_sequence || par.cap ? &ctthreshold : NULL
         ;  dim tblidx = 0
                                                      /* 4 extra because of ctwindow.n_***
                                                       * hack further below.
                                                      */
         ;  if (!(dt.tablect = mycalloc(n_analyse * (dt.wordlist.n_words+4), sizeof dt.tablect[0])))
            return mainfail(&dt, "could not alloc table")

         ;  if (acc_init(&ctwindow, par.K_g, dt.hash_template_g))
            return mainfail(&dt, NULL)

         ;  if (INTERFACE(par, MODE_GAPPY) && acc_init(&ctcache, par.K_g, dt.hash_template_g))
            return mainfail(&dt, NULL)

         ;  if (ctacc_threshold && acc_init(ctacc_threshold, par.K_g, dt.hash_template_g))
            return mainfail(&dt, NULL)

         ;  for (i=0;i<n_analyse;i++)
            {  unit* ut = dt.thecollection->list+i
            ;  if (par.lastbin && i >= par.lastbin)
               break

            ;  acc_reset(&ctwindow)

            ;  status
               =  INTERFACE(par, MODE_GAPPY)
               ?  acc_addto_gapped(par.K_g, ut, &ctwindow, par.gap_min, par.gap_max, par.gap_gap, &ctcache)
               :  acc_addto(par.K_g, ut, &ctwindow, INTERFACE(par, MODE_R2ON))
            ;  if (status)
               return mainfail(&dt, "addto count")

            ;  acc_iter_reset(&iterator)

            ;  while (acc_iterate(&iterator, &ctwindow, &status)) 
               dt.tablect[tblidx++]  =  iterator.count1

            ;  if (status)
               return mainfail(&dt, "error while iterating over counts")

            ;  dt.tablect[tblidx++]   =  ctwindow.n_bases
            ;  dt.tablect[tblidx++]   =  ctwindow.n_mask
            ;  dt.tablect[tblidx++]   =  ctwindow.n_syl
            ;  dt.tablect[tblidx++]   =  ctwindow.n_unsyl
         ;  }

            {  char sylbuf[16]
            ;  dim n_bin_used = par.lastbin && par.lastbin < n_analyse ? par.lastbin : n_analyse
            ;  dim n_passed = 0
            ;  if (MAIN_MODE(par, MAIN_FILTER))
               {  for (i=0;i<n_bin_used;i++)
                  {  dim j, sum = 0
                  ;  int ok = 1
                  ;  for (j=0;j<dt.wordlist.n_words;j++)
                     {  dim ct = dt.tablect[i*(dt.wordlist.n_words+4) + j]
                     ;  if (INTERFACE(par, MODE_SUM))
                        sum += ct
                     ;  else
                        {  ok =
                           (  (par.unit_at_most  < 0 || ct <= par.unit_at_most )
                           && (par.unit_at_least < 0 || ct >= par.unit_at_least)
                           )  ?
                           1  :  0
                        ;  if ( (!ok && INTERFACE(par, MODE_AND)) || (ok && !INTERFACE(par, MODE_AND)) )
                           break
                     ;  }
                     }
                     if (INTERFACE(par, MODE_SUM))
                     ok =
                     (  (par.unit_at_most  < 0 || sum <= par.unit_at_most )
                     && (par.unit_at_least < 0 || sum >= par.unit_at_least)
                     )  ?
                     1  :  0
                  ;  if (ok)
                        fputs(dt.thecollection->list[i].sid, par.fpresult)
                     ,  fputc('\n', par.fpresult)
                     ,  n_passed++
               ;  }
                  tell(VERBOSE_ONCE, "sequences passed %lu", (unsigned long) n_passed)
            ;  }
               else if (MAIN_MODE(par, MAIN_COUNT))
               {  fputs("id", par.fpresult)
               ;  for (i=0;i<dt.wordlist.n_words;i++)
                     fputc('\t', par.fpresult)
                  ,  fputs(get_sylmer(par.K_g, dt.wordlist.words ? dt.wordlist.words[i] : i, sylbuf), par.fpresult)
               ;  fputs("\tbases\tmasked\tsyl\tunsyl\trank\n", par.fpresult)
               ;  for (i=0;i<n_bin_used;i++)
                  {  dim j
                  ;  fputs(dt.thecollection->list[i].sid, par.fpresult)
                  ;  for (j=0;j<dt.wordlist.n_words+4;j++)
                     fprintf(par.fpresult, "\t%lu", dt.tablect[i*(dt.wordlist.n_words+4) + j])
                  ;  fprintf(par.fpresult, "\t%lu", (unsigned long) (i+1))
                  ;  fputc('\n', par.fpresult)
               ;  }
               }
            }
                  /* fixme; these free/close instances are not yet designed to follow error paths
                   * solution: put them together in a new data bundle and introduce
                   * a callback routine in dataframe struct.
                  */
            myfclose(&par.fpresult)
         ;  acc_free(&ctwindow)
         ;  if (INTERFACE(par, MODE_GAPPY))
            acc_free(&ctcache)
         ;  if (ctacc_threshold)
            acc_free(ctacc_threshold)
         ;  return myfinish(&dt, 0)
      ;  }
#endif

