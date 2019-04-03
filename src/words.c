

#include "inttypes.h"
#include "words.h"
#include "sylio.h"
#include "acc.h"
#include "util.h"
#include "interface.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>


/* TODO:

 * sanity check that wl->n_words is kept consistent at all times.
 ------> provide setter routine, or make wlist_extend do it.

 * the dyad routines assumes wl is empty.
*/


sylstatus wlist_malloc(wlist* wl, dim n_amount)
   {  unsigned* new
   ;  if (n_amount <= wl->n_words_alloc)
      return SYL_SUCCESS
   ;  if (!(new = myrealloc(wl->words, n_amount * sizeof wl->words[0])))
      return SYL_FAIL
   ;  wl->words = new
   ;  wl->n_words_alloc = n_amount
   ;  return SYL_SUCCESS
;  }


sylstatus wlist_extend(wlist* wl)
   {  dim n_new = wl->n_words_alloc ? wl->n_words_alloc * 1.41 : 4096
   ;  if (wl->n_words < wl->n_words_alloc)
      return SYL_SUCCESS
   ;  return wlist_malloc(wl, n_new)
;  }


void wlist_free(wlist* wl)
   {  if (wl->words)
      free(wl->words)
   ;  wl->words = NULL
   ;  wl->n_words = 0
   ;  wl->n_words_alloc = 0
;  }


sylstatus wlist_full(wlist* wl, unsigned k)
   {  dim i
   ;  sylstatus status = wlist_malloc(wl, LSIZE(k))
   ;  if (status)
      return status
   ;  for (i=0;i<LSIZE(k);i++)
      wl->words[i] = i   
   ;  wl->n_words = LSIZE(k)
   ;  return SYL_SUCCESS
;  }


static sylstatus expand_iterate(unsigned char* o, unsigned char* g, int p, int l, wlist* wl)
   {  int i
   ;  if (p == l)
      {  unsigned kmer = 0
      ;  if (wlist_extend(wl))
         return SYL_FAIL
      ;  for (i=0;i<l;i++)
         {  int c = o[i]
         ;  kmer <<= 2
                        /* bit 1 -> A { 0 }, 2 -> C { 1 }, 4 -> G { 2 }, 8 -> T { 3 } */
         ;  kmer |= (c == 1 ? 0 : c == 2 ? 1 : c == 4 ? 2 : c == 8 ? 3 : 0)
      ;  }
         wl->words[wl->n_words] = kmer
      ;  wl->n_words++
   ;  }
      else
      {  unsigned c = g[p]
      ;  for (i=0;i<4;i++)
         {  if ((1 << i) & c)
            {  o[p] = 1 << i
            ;  if (expand_iterate(o, g, p+1, l, wl))
               return SYL_FAIL
         ;  }
         }
      }
      return SYL_SUCCESS
;  }


static dim generator_fill
(  unsigned char* generator
,  dim g_len
,  unsigned K_ref
,  const char* pattern
)
   {  dim i = 0, g = 0
   ;  dim l = strlen(pattern)
   ;  int in_range = 0
   ;  const char* error = NULL
   ;  while (i < l)
      {  if (g >= g_len)
            fprintf(stderr, "expression too large\n")
         ,  exit(1)
      ;  if (in_range)
         {  switch(pattern[i])
            {  case ']' :   in_range = 0; g++;  break    /* noteme: coulds use |= 1 << basemap[pattern[i]] in all cases */
            ;  case 'A' :   generator[g] |= 1;  break
            ;  case 'C' :   generator[g] |= 2;  break
            ;  case 'G' :   generator[g] |= 4;  break
            ;  case 'T' :   generator[g] |= 8;  break
            ;  default  :   error  = "unexpected range character"
         ;  }
            i++
         ;  if (error)
            break
         ;  continue
      ;  }

         switch(pattern[i])
         {  case '[' :     in_range = 1; generator[g] = 0;  break
         ;  case '.' :     generator[g] = 0xF;  break      /* four bits set */

         ;  case 'A' :     generator[g] = 1;    break
         ;  case 'T' :     generator[g] = 8;    break
         ;  case 'U' :     generator[g] = 8;    break
         ;  case 'W' :     generator[g] = 9;    break      /* Weak, A or T  */

         ;  case 'C' :     generator[g] = 2;    break
         ;  case 'G' :     generator[g] = 4;    break
         ;  case 'S' :     generator[g] = 6;    break      /* Strong, C or G */

         ;  case 'K' :     generator[g] = 12;   break      /* Keto, G or T   */
         ;  case 'M' :     generator[g] = 3;    break      /* aMino, A or C  */
         ;  case 'R' :     generator[g] = 5;    break      /* puRine, A or G */
         ;  case 'Y' :     generator[g] = 10;   break      /* pYrimidine, C or T */

         ;  case 'B' :     generator[g] = 14;   break      /* not A           */
         ;  case 'V' :     generator[g] = 7;    break      /* not T/U         */
         ;  case 'H' :     generator[g] = 11;   break      /* not G           */
         ;  case 'D' :     generator[g] = 13;   break      /* not C           */

         ;  case '*' :     while(g > 0 && g < K_ref)
                           {  generator[g] = generator[g-1]
                           ;  g++
                           ; 
                           }
                           break
         ;  default  :   error  = "unexpected pattern character"
      ;  }
         if (!in_range && pattern[i] != '*')    /* fixme: '*' was hacked in */
         g++
      ;  if (error)
         break
      ;  i++
   ;  }

      if (error)
      {  miaow("%s", error)
      ;  return 0
   ;  }
      return g
;  }


static sylstatus parse_pattern(const char* pattern, unsigned* K_found, unsigned K_ref, wlist* wl)
   {  unsigned char generator[15]
   ;  unsigned char output[15]
   ;  int g

;if(0)fprintf(stderr, "pat [%s] K_ref %d\n", pattern, (int) K_ref)

   ;  if (!(g = generator_fill(generator, 15, K_ref, pattern)))
      return SYL_FAIL

   ;  if (K_ref && g != K_ref)
      {  miaow("pattern expand -- conflicting K value (expansion: %d)", (int) g)
      ;  return SYL_FAIL
   ;  }

      memcpy(output, generator, g)
   ;  if (expand_iterate(output, generator, 0, g, wl))
      wlist_free(wl)

;if(0)fprintf(stderr, "%d words parsed first %d\n", (int) wl->n_words, (int) wl->words[0])
   ;  *K_found = g
   ;  return SYL_SUCCESS
;  }


static sylstatus parse_words(const char* list, wlist* wl, unsigned* K_found, unsigned K_ref)
   {  const char* s = list, *se = list+strlen(list)
   ;  char pattern[128]
   ;  sylstatus status

   ;  if ((status = wlist_extend(wl)))
      return status

   ;  do
      {  dim l
      ;  s = strchr(list, ',')
      ;  l = (s ? s : se) -list
      ;  status = SYL_FAIL

      ;  if (l > 127)
         {  miaow("pattern too large")
         ;  break
      ;  }
         if (!l)
         {  miaow("emptiness has a bound")
         ;  break
      ;  }

         memcpy(pattern, list, l)
      ;  pattern[l] = '\0'

      ;  if (parse_pattern(pattern, K_found, K_ref, wl))
         break

      ;  if (!K_found[0])
         {  miaow("emptiness knows no bound")
         ;  break
      ;  }

         if (!K_ref)
         K_ref = K_found[0]
      ;  list = s+1
      ;  status = SYL_SUCCESS
   ;  }
      while (s)

   ;  return status
;  }


sylstatus get_bypattern(const char* list, wlist* wl, unsigned* K_found)
   {  unsigned K_ref = K_found[0]
   ;  return parse_words(list, wl, K_found, K_ref)
;  }


static sylstatus read_words (INPUT fp, wlist* wl, unsigned* K, float** frequencies)
   {  char* line
   ;  dim linelen       =  0
   ;  dim linecount     =  0
   ;  dim n_words       =  0
   ;  float* frequencies_realloc
   ;  unsigned mask     =  0
   ;  unsigned k        =  0
   ;  sylstatus status  =  SYL_SUCCESS
   ;  dim frequencies_n_alloc  =  0

   ;  if ((status = wlist_extend(wl)))
      return status

   ;  if (frequencies && frequencies_n_alloc < wl->n_words_alloc)
      {  frequencies_realloc
         =  *frequencies
         =  mymalloc(wl->n_words_alloc * sizeof frequencies[0])
      ;  frequencies_n_alloc = wl->n_words_alloc
   ;  }

      if (frequencies && !*frequencies)
      return SYL_FAIL

   ;  while ((line = fread_line(fp, NULL, 0, &linelen, &linecount, &status)) && !status)
      {  size_t l = strcspn(line, " \t\r")
      ;  dim j

      ;  status = SYL_FAIL       /* success if we reach end of block */
      ;  line[l] = '\0'

      ;  if (!n_words)
         k = l
      ;  else if (l != k)
         {  miaow("inconsistent word length at line %d", (int) linecount)
         ;  break
      ;  }
         if ((status = wlist_extend(wl)))
         break

      ;  if (frequencies && frequencies_n_alloc < wl->n_words_alloc)
         {  frequencies_realloc
         =  myrealloc(*frequencies, wl->n_words_alloc * sizeof frequencies[0])
         ;  if (!frequencies_realloc)
            break
         ;  *frequencies = frequencies_realloc
         ;  frequencies_n_alloc = wl->n_words_alloc
      ;  }

         for (j=0;j<l;j++)
         line[j] = toupper((unsigned char) line[j])

      ;  wl->words[wl->n_words] = sylid_from_buf(line, k, &mask)

      ;  if (mask)
         {  miaow("funny word %s at line %d", line, (int) linecount)
         ;  break
      ;  }

         if (frequencies)
         {  float f
         ;  if
            (  l == linelen
            || 1 != sscanf(line+l+1, "%f", &f)
            || f < 0 || f > 1
            )
            {  miaow("expect frequency at line %d", (int) linecount)
            ;  break
         ;  }
            (*frequencies)[wl->n_words] = f
      ;  }
         wl->n_words++
      ;  free(line)
      ;  status = SYL_SUCCESS
   ;  }

      if (status == SYL_DONE)
      status = SYL_SUCCESS
   ;  else
      {  if (frequencies && *frequencies)
         free(*frequencies)
      ;  if (line)
         free(line)
      ;  status = SYL_FAIL
   ;  }
      *K = k
   ;  return status
;  }


int cmp_unsigned(const void* p1, const void* p2)
   {  unsigned u1 = * (unsigned*) p1, u2 = * (unsigned*) p2
   ;  return u1 > u2 ? 1 : u1 < u2 ? -1 : 0
;  }



void wlist_remove
(  wlist*   wl
,  wlist*   wldelete
)
   {  dim  offset, n_found = 0
   ;  qsort(wldelete->words, wldelete->n_words, sizeof wldelete->words[0], cmp_unsigned)
   ;  for (offset=0;offset<wl->n_words;offset++)
      {  if (bsearch(wl->words+offset, wldelete->words, wldelete->n_words, sizeof wldelete->words[0], cmp_unsigned))
         n_found++
      ;  else if (n_found)
         wl->words[offset-n_found] = wl->words[offset]
   ;  }
      wl->n_words -= n_found
;  }


dim wlist_uniq(wlist* wl)
   {  dim offset, n_dup = 0
   ;  qsort(wl->words, wl->n_words, sizeof wl->words[0], cmp_unsigned)
   ;  for (offset=1;offset<wl->n_words;offset++)
      {  if (wl->words[offset] == wl->words[offset-n_dup-1])
         n_dup++
      ;  else if (n_dup)
         wl->words[offset-n_dup] = wl->words[offset]
   ;  }
      wl->n_words -= n_dup
   ;  return n_dup
;  }


int cmp_uchar(const void* p1, const void* p2)
   {  const unsigned char u1 =  *(unsigned char*) p1, u2 = *(unsigned char*) p2
   ;  return basemap[u1] < basemap[u2] ? -1 : basemap[u1] > basemap[u2] ? 1 : 0
;  }


static sylstatus empty_bag(char* kmer, unsigned k, unsigned *bag, harray* ha, wlist* wl, dim n_left)
   {  unsigned mask = 0
   ;  if (!n_left)
      {  unsigned w = sylid_from_buf(kmer, k, &mask)
      ;  int found = 0
      ;  HISPRESENT(found, ha, w)
      ;  if (!found)
         {  if (wlist_extend(wl))
            return SYL_FAIL
         ;  wl->words[wl->n_words++] = w
         ;  HSET(ha, w, 1)
      ;  }
/* else fprintf(stderr, "found %s\n", kmer); */
      }
      else
      {  dim p = 0, j = 0
      ;  while (kmer[p] != 0)
         p++
      ;  for (j=0;j<4;j++)
         {  if (bag[j] > 0)
            {  kmer[p] = j == 0 ? 'A' : j == 1 ? 'C' : j == 2 ? 'G' : 'T'
            ;  bag[j]--
            ;  empty_bag(kmer, k, bag, ha, wl, n_left-1)
            ;  bag[j]++
            ;  kmer[p] = 0
         ;  }
         }
      }
      return SYL_SUCCESS
;  }


sylstatus get_bags(wlist* wl, unsigned k)
   {  char kmer[16] = { 0 }
   ;  harray* ha
   ;  dim i, n_words_old
   ;  unsigned mask = 0

   ;  wlist_uniq(wl)
   ;  n_words_old = wl->n_words
                              /* first add bags for all the current words */
   ;  for (i=0;i<n_words_old;i++)
      {  get_sylmer(k, wl->words[i], kmer)
      ;  qsort(kmer, k, 1, cmp_uchar)  /* create canonical (sorted) bag for this word */
      ;  if (wlist_extend(wl))
         return SYL_FAIL
      ;  wl->words[wl->n_words++] = sylid_from_buf(kmer, k, &mask)
      ;  if (mask)
         {  miaow("unexpected mask while bagging")
         ;  return SYL_FAIL
      ;  }
      }
      tell
      (  VERBOSE_ONCE
      ,  "[bag] added bags for %d unique words"
      ,  (int) n_words_old
      )

                              /* move all the bags to front */
   ;  memmove(wl->words, wl->words+n_words_old, sizeof wl->words[0] * (wl->n_words - n_words_old))
                              /* set the length correctly   */
   ;  wl->n_words -= n_words_old
                              /* then make the list unique */
   ;  n_words_old = wl->n_words
   ;  wlist_uniq(wl)
   ;  tell
      (  VERBOSE_ONCE
      ,  "[bag] removed %d redundant bags list now has %d"
      ,  (int) (n_words_old - wl->n_words)
      ,  (int) wl->n_words
      )
                              /* and hash all the words in it */
   ;  ha = harray_new(wl, LOMEGA(k), k)
                              /* make it blow up if we use words member,
                               * as we change it underneath. (this is not very violent,
                               * this member should serve only specific iterators).
                              */
   ;  ha->words = NULL
   ;  n_words_old = wl->n_words

   ;  for (i=0;i<n_words_old;i++)
      {  unsigned bag[4] = { 0 }
      ;  dim j
      ;  get_sylmer(k, wl->words[i], kmer)         /* write word in buffer */
      ;  for (j=0;j<k;j++)
         bag[basemap[(unsigned char) kmer[j]]]++   /* create counts for all the bases */
      ;  memset(kmer, 0, k)                        /* make buffer available for bags */
      ;  if (empty_bag(kmer, k, bag, ha, wl, k))   /* recurse */
         return SYL_FAIL
   ;  }
      harray_free(ha)
   ;  tell
      (  VERBOSE_ONCE
      ,  "[bag] inserted %d new words"
      ,  (int) (wl->n_words - n_words_old)
      )

   ;  wlist_uniq(wl)          /* sort it, (should already be unique) */
   ;  return SYL_SUCCESS
;  }


dim wlist_entropy(wlist* wl, unsigned K, double e)
   {  dim offset, n_rm = 0
   ;  for (offset=0;offset<wl->n_words;offset++)
      {  if (get_word_entropy(wl->words[offset], K) < e)
         n_rm++
      ;  else if (n_rm)
         wl->words[offset-n_rm] = wl->words[offset]
   ;  }
      wl->n_words -= n_rm
   ;  return n_rm
;  }


            /* the error checking in this code is not the summit of elegance */
sylstatus get_wordlist
(  INPUT fpwords           /* closed by caller */
,  wlist* wl
,  float** oracle_pp
,  unsigned* K_pp          /* check for consistency if *K_pp != 0, otherwise set it */
)
   {  unsigned K_found = 0
   ;  unsigned K_ref = *K_pp
   ;  sylstatus status = SYL_SUCCESS

   ;  do
      {  status = read_words(fpwords, wl, &K_found, oracle_pp)
      ;  if (status)
         {  miaow("[get_wordlist] could not read wordlist from file")
         ;  break
      ;  }
         if (K_ref && K_found != K_ref)
         {  miaow("[get_wordlist] -k value conflicts with -word contents")
         ;  status = SYL_FAIL
         ;  break
      ;  }
         *K_pp = K_found
   ;  }
      while (0)
   ;  return status
;  }


sylstatus get_dyadlist
(  unsigned modes
,  wlist* wl
,  unsigned Kh             /* dyad length (half of K). */
)
   {  int repeat = (modes & DYAD_REPEAT) ? 1 : 0
   ;  int invert = (modes & DYAD_INVERT) ? 1 : 0
   ;  unsigned* words = mymalloc((invert + repeat) * LSIZE(Kh) * sizeof words[0])
   ;  unsigned i_repeat, i_invert = 0

   ;  if (wlist_malloc(wl, (invert + repeat) * LSIZE(Kh)))
      return SYL_FAIL

   ;  for (i_repeat=0;i_repeat<LSIZE(Kh);i_repeat++)
      {  if (repeat)
         wl->words[i_repeat] = i_repeat << (Kh * 2) | i_repeat
      ;  if (invert)
         {  dim j
         ;  unsigned rc = 0
         ;  for (j=0;j<Kh;j++)
            rc |= (~(i_repeat >> 2*j) & 3) << 2*(Kh-j-1)
         ;  if (rc != i_repeat)
            wl->words[repeat * LSIZE(Kh) + i_invert++] = (i_repeat << (Kh * 2)) | rc
      ;  }
      }
      wl->n_words = repeat * LSIZE(Kh) + i_invert
   ;  return SYL_SUCCESS
;  }



double get_word_entropy(unsigned w, unsigned K)
   {  unsigned ct[4]
   ;  double e = 0.0
   ;  dim i
   ;  if (!K)
      return 0.0
   ;  memset(ct, 0, 4 * sizeof ct[0])
   ;  for (i=0;i<K;i++)
      ct[(w >> (2*i)) & 3]++
   ;  for (i=0;i<4;i++)
      {  double p = ct[i] * 1.0 / K
      ;  if (p)
         e -= p * log(p) / log(4)
   ;  }
      return e
;  }


double get_string_entropy (const char* word)
   {  unsigned ct[4]
   ;  dim i, l = strlen(word)
   ;  double e = 0.0
   ;  if (!l)
      return 0.0
   ;  memset(ct, 0, 4 * sizeof ct[0])
   ;  for (i=0;i<l;i++)
      ct[basemap[(unsigned char) word[i]]]++
   ;  for (i=0;i<4;i++)
      {  double p = ct[i] * 1.0 / l
      ;  if (p)
         e -= p * log(p) / log(4)
   ;  }
      return e
;  }




