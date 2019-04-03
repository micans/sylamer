

#include "unit.h"
#include "util.h"
#include "query.h"
#include "acc.h"
#include "interface.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

struct mask_stats
{  unsigned n_nucs[16]
;  unsigned n_seqs[16]
;  unsigned n_pats[16]
;
}  ;

static const char* fnrepeat_g = NULL;
static FILE* fprepeat_g = NULL;


int unit_cmp_by_length(const void* p1, const void* p2)
   {  const unit *ut1 = p1, *ut2 = p2;
   ;  return ut1->n_bases < ut2->n_bases ? -1 : ut1->n_bases == ut2->n_bases ? 0 : 1
;  }


int unit_cmp_by_revlength(const void* p1, const void* p2)
   {  const unit *ut1 = p1, *ut2 = p2;
   ;  return ut1->n_bases > ut2->n_bases ? -1 : ut1->n_bases == ut2->n_bases ? 0 : 1
;  }


int unit_cmp_by_rank(const void* p1, const void* p2)
   {  const unit *ut1 = p1, *ut2 = p2;
   ;  return ut1->rank < ut2->rank ? -1 : ut1->rank == ut2->rank ? 0 : 1
;  }


int unit_cmp_by_revrank(const void* p1, const void* p2)
   {  const unit *ut1 = p1, *ut2 = p2;
   ;  return ut1->rank > ut2->rank ? -1 : ut1->rank == ut2->rank ? 0 : 1
;  }


static void unit_free(unit* ut)
   {  if (ut->bases)    free(ut->bases)
   ;  if (ut->sid)      free(ut->sid)
;  }


#if 0
static void unit_init(unit* ut)
   {  ut->bases   =  NULL
   ;  ut->annot   =  ""
   ;  ut->n_bases =  0
   ;  ut->sid     =  NULL
   ;  ut->rank    =  0
;  }
#endif


dim sequence_normalize
(  unit* ut
)
   {  unsigned char* bases = ut->bases
   ;  dim  n_bases = ut->n_bases, n_masked = 0
   ;  dim j, n_weird = 0
   ;  double n_gc = 0.0

   ;  for (j=0;j<n_bases;j++)
      {  bases[j] = toupper(bases[j])
      ;  if (basemap[bases[j]] == 4)
         n_masked++
      ;  else if (basemap[bases[j]] > 4)
         n_weird++
      ;  else if (basemap[bases[j]] == 1 || basemap[bases[j]] == 2)
         n_gc++
   ;  }
      ut->n_masked = n_masked
   ;  ut->gccontent = (n_bases - n_masked -n_weird) ? n_gc * 1.0 / (n_bases - n_masked - n_weird) : 0.5
   ;  return n_weird
;  }


void collection_init
(  collection* cl
)
   {  cl->list = NULL
   ;  cl->list_n = 0
;  }


void collection_release(collection* cl)
   {  dim i
   ;  if (cl->list)
      {  for (i=0;i<cl->list_n;i++)
         unit_free(cl->list+i)
      ;  free(cl->list)
   ;  }
      cl->list_n = 0
   ;  cl->list = NULL
;  }

#define REVERSE_COMPLEMENT(b) \
 b == 'A' ? 'T' : b == 'C' ? 'G' :  b == 'G' ? 'C' : b == 'T' ? 'A' : b  /* e.g. 'N' or 'X' */



   /* this leaks non-trivial memory under memory shortage,
    * because in the construction in the main loop new elements
    * are inserted from the (new) end and cl->list_n is not
    * capable of tracking this.
    * Fixable with some thought. todo/fixme.
   */
sylstatus collection_revcompl(collection* cl)
   {  dim i
   ;  unit* list_realloc = myrealloc(cl->list, 2 * cl->list_n * sizeof cl->list[0])

   ;  if (!list_realloc)
      return SYL_FAIL

   ;  cl->list = list_realloc

   ;  for (i=0;i<cl->list_n;i++)
      {  dim j    =  cl->list_n-i-1
      ;  unit* src=  cl->list+j
      ;  unit* dst=  cl->list+j+j
      ;  unsigned k
                  /* we work from last to first.
                   * each unit yields two new units, stored consecutively.
                  */

      ;  dst[0] =  src[0]
      ;  dst[0].rank = 2 * src[0].rank

      ;  dst[1].bases = mymalloc((dst[0].n_bases+1) * sizeof dst[1].bases[0])
      ;  dst[1].sid = strdup(dst[0].sid)

      ;  if (!dst[1].bases || !dst[1].sid)
         break

      ;  dst[1].n_bases = dst[0].n_bases
      ;  dst[1].n_masked = dst[0].n_masked
      ;  dst[1].gccontent = dst[0].gccontent
      ;  dst[1].rank = dst[0].rank+1

      ;  dst[1].bases[dst[1].n_bases] = '\0'

      ;  for (k=0;k<dst[1].n_bases;k++)
         {  int b = dst[0].bases[dst[0].n_bases-k-1]
         ;  dst[1].bases[k] = REVERSE_COMPLEMENT(b)
      ;  }
   ;  }

      if (i != cl->list_n)       /* fixme; try to get cl->list_n correct */
      return SYL_FAIL

   ;  cl->list_n *= 2
   ;  return SYL_SUCCESS
;  }



sylstatus collection_revcompl2(collection* cl, dim gap_max)
   {  dim i

   ;  for (i=0;i<cl->list_n;i++)
      {  unit* ut =  cl->list+i
      ;  dim N = ut->n_bases, N_new = 2 * N + 1 + gap_max
      ;  dim k
                              /* + 1 for '\0', gap_max + 1 for separating N */
      ;  unsigned char* dbase = mymalloc((N_new + 1) * sizeof ut->bases[0]);

      ;  if (!dbase)
         return SYL_FAIL         /* note that cl is still in consistent state */

      ;  memcpy(dbase, ut->bases, N)
      ;  memset(dbase+N, 'N', gap_max+1)
      ;  dbase[N_new] = '\0'

      ;  for (k=0;k<N;k++)
         {  int b = ut->bases[N-k-1]
         ;  dbase[N+k+1+gap_max] = REVERSE_COMPLEMENT(b)
      ;  }
         free(ut->bases)
      ;  ut->bases = dbase
      ;  ut->n_bases = N_new
      ;  ut->n_masked *= 2
   ;  }
      return SYL_SUCCESS
;  }


void collection_shuffle(collection* cl, dim n_universe)
   {  dim n = n_universe
   ;  while (n > 0)
      {  unsigned long r = (random() >> 3) % n             /* Fisher-Yates shuffle */
      ;  unit src  =  cl->list[n-1]
      ;  unit dst  =  cl->list[r]
      ;  cl->list[r]  =  src
      ;  cl->list[n-1]=  dst
      ;  n--
   ;  }
   }


   /* keep the remainder of the repetitive stretch, remove cycles in front. */
static int splice_regular
(  unsigned char* bases
,  int   o
,  int   ival_len
,  int   cyclen
,  int   K
)  
   {  int keep = K + cyclen - 1, i, rm = 0

   ;  rm = cyclen * ((ival_len - keep) / cyclen)

;if(0)fprintf(stderr, "bases offset %d length %d rm %d [%s]\n",
   (int) o, (int) ival_len, (int) rm, bases)

   ;  for (i=0;i<rm;i++)
      bases[o+i] = '@'

;if(0)fprintf(stderr, "bases offset %d length %d rm %d [%s]\n",
   (int) o, (int) ival_len, (int) rm, bases)
   ;  return rm
;  }


   /* the repetitive bit is spliced in two, with 'keep' bases remaining
    * at each end, separated by a stretch of Ns which is a multiple of cyclen.
   */
static int splice_gap
(  unit* ut
,  int   o
,  int   ival_len
,  int   cyclen
,  int   K
,  int   gap_size
,  int   * skip
)
   {  int keep = K + cyclen -1, i
   ;  int delta = 2 * keep + gap_size - ival_len
   ;  int rm = ival_len - 2 * keep

   ;  if (rm < 0)
      rm = 0

   ;  if (delta > 0)
      {  if (!(ut->bases = realloc(ut->bases, ut->n_bases + delta + 1)))
            miaow("panic memory failure")
         ,  exit(1)           /* TODO fixme escalate to sylmain() */
;if(0)fprintf(stdout, "delta %d keep %d\n", delta, keep)
      ;  memmove
         (  ut->bases+o+ival_len+delta-keep
         ,  ut->bases+o+ival_len-keep
         ,  ut->n_bases - o - ival_len + keep
         )
      ;  ut->n_bases += delta         
      ;  ut->bases[ut->n_bases] = '\0'
      ;  ival_len += delta
   ;  }
;if(0)fprintf(stderr, "rm %d keep %d ival_len %d\n", (int) rm, (int) keep, (int) ival_len);
      for (i=keep;i<ival_len-keep;i++)
      ut->bases[o+i] = rmchar_G
   ;  skip[0] = ival_len -keep + 1
   ;  return rm
;  }


dim sequence_contract
(  unsigned char* bases
,  dim n_bases
,  unsigned c
)
   {  dim i = 0, d = 0
   ;  while (i<n_bases)
      {  if (bases[i] != c)
         {  if (d < i)
            bases[d] = bases[i]
         ;  d++
      ;  }
         i++
   ;  }
      bases[d] = '\0'
   ;  return d
;  }



static dim sequence_rm
(  unit* ut
,  int max_cycle_length
,  int gap_size_required
,  int theK
,  unsigned maskall
,  struct mask_stats *ms
)
   {  dim i = 0
   ;  unsigned char* bases = ut->bases
   ;  dim n_bases = ut->n_bases
   ;  dim rm = 0
   ;  unsigned seq_count[16] = { 0 }

   ;  while(i < n_bases)
      {  int c = 0, c_stretch = 0, r
      ;  dim j
      ;  for (r=1; r<=max_cycle_length;r++)
         {  j = i+r
         ;  while
            (  j<n_bases
            && basemap[(unsigned char) bases[j]] < 4
            && (unsigned char) bases[j] == bases[j-r]
            )
            j++
         ;  if (j-i > c_stretch)
               c_stretch = j-i
            ,  c = r
      ;  }
         if (c_stretch >= theK + 2 * c)
         {  int skip = 0, rmed = 0
;if(0)fprintf(stderr, "hit %d %.*s\n", (int) i, (int) c, bases+i)

         ;  if (fprepeat_g)
            fprintf
            (  fprepeat_g
            ,  ">%d-%d-%s\n%.*s\n"
            ,  (int) c
            ,  (int) (seq_count[c] + 1)
            ,  ut->sid
            ,  c_stretch
            ,  bases+i
            )

         ;  if (maskall)
            {  memset(bases+i, rmchar_G, c_stretch)
            ;  skip = c_stretch
            ;  rm += (rmed = c_stretch)
         ;  }
            else if (gap_size_required)
            {  rm += (rmed = splice_gap(ut, i, c_stretch, c, theK, gap_size_required, &skip))
            ;  bases = ut->bases    /* might have been realloc'ed */
            ;  n_bases = ut->n_bases
         ;  }
            else                    /* splice_reg inserts _ for sequence that should be deleted */
            {  rm += (rmed = splice_regular(bases, i, c_stretch, c, theK))
            ;  skip = c_stretch - theK + 1
         ;  }
;if(0)fprintf(stderr, " cycle %d length %d skip %d rm %d %.6s\n", (int) c, (int) c_stretch, (int) rm, (int) skip, bases+i+skip)
         ;  i += skip
         ;  ms->n_nucs[c] += rmed
         ;  if (rmed)
            seq_count[c]++
      ;  }
         else
         i++
   ;  }
      if (!gap_size_required)
         ut->n_bases = sequence_contract(bases, n_bases, '@')
      ,  ut->bases[ut->n_bases] = '\0'
   ;  for (i=0;i<16;i++)
         ms->n_seqs[i] += seq_count[i] ? 1 : 0
      ,  ms->n_pats[i] += seq_count[i]
   ;  return rm
;  }


static dim validate_sequence
(  unit* ut
,  int warn_funny
,  int length_co
,  int max_cycle_length
,  int gap_size_required
,  int theK
,  unsigned modes
,  struct mask_stats *ms
)
   {  dim n_weird = sequence_normalize(ut)
   ;  dim lco = length_co > 0 ? length_co : -length_co
   ;  dim rm = 0
   ;  unsigned maskall = modes & MODE_MASKALL
   ;  unsigned ulterior = modes & MODE_ULTERIOR

   ;  if (warn_funny && n_weird)
      miaow
      (  "warning: found %lu unrecognized base symbols for <%s>"
      ,  (unsigned long) n_weird
      ,  ut->sid
      )
   ;  if (lco && ut->n_bases > lco)
      {  if (ulterior && 2 * lco < ut->n_bases)
         memset(ut->bases, 'N', ut->n_bases - 2 * lco) 
      ;  else
         {  if (length_co < 0)
            memmove(ut->bases, ut->bases+ut->n_bases - lco, lco)
         ;  ut->n_bases = lco
         ;  ut->bases[lco] = '\0'
      ;  }
      }
      if (max_cycle_length)
      rm = sequence_rm(ut, max_cycle_length, gap_size_required, theK, maskall, ms)
   ;  return rm
;  }


static int search_id
(  char* sid
,  query* thequery
,  dim* rank_found
,  dim  linecount
)
   {  qel* match =  NULL
   ;  *rank_found = -1u

                           /* search for sid in query file */
   ;  if (thequery->qel_list)
      {  qel dummy
      ;  dummy.sid = sid
      ;  match
         =  thequery->qel_n_uniq
         ?  bsearch       /* we should probably move this to query.c */
            (  &dummy
            ,  thequery->qel_list
            ,  thequery->qel_n_uniq
            ,  sizeof thequery->qel_list[0]
            ,  qel_cmp_by_id
            )
         :  NULL           /* bsearch usually crashes on empty list :-( */

                           /* if found, set rank to query rank */
      ;  if (match && !match->found)
         {  match->found++
         ;  *rank_found = match->rank
         ;  thequery->qel_n_found++
      ;  }
         else if (match && match->found)
         {  miaow
            (  "warning: ignoring duplicate id %s (query rank %d) in fasta file line %u"
            ,  match->sid
            ,  (int) match->rank
            ,  (unsigned) linecount
            )
         ;  match = NULL
      ;  }
      }
      return match ? 1 : 0
;  }
                              /* fixme: if no query, we do not
                               * detect duplicates in the fasta file
                               * at the moment - as with non-found duplicates.
                              */



static char* emptyseq(void)
   {  char* seq = mymalloc(1)
   ;  if (seq)
      seq[0] = '\0'
   ;  return seq
;  }


   /* Created this routine because valgrind (on optimized code compile without -g)
    * gave errors/warnings during a qsort complaining about unitialized values.
    * Checking with collection_cat gave no such errors, so I assume it
    * was a valgrind artefact.
   */
void collection_cat(collection* cl)
   {  dim i
   ;  for (i=0;i<cl->list_n;i++)
      {  unit* ut = cl->list+i
      ;  fprintf(stdout, "%s [%s] >%s< %d %d\n", ut->sid, ut->annot, ut->bases, (int) ut->n_bases, (int) ut->rank)
   ;  }
   }


struct read_state
{  char*    cache_line
;  dim      cache_length
;  dim      linecount
;  dim      n_skip
;  dim      rank_notfound
;
}  ;


   /* search ID.
    * decide: flush or not?
    * if flush  -> flush, return IGNORE or DONE or FAIL
    * if !flush -> set sid, annot, rank, return SUCCESS
   */
sylstatus read_id
(  INPUT       fp
,  unit*       dest
,  query  *    thequery
,  int         universe
,  struct read_state * state
)
   {  size_t l = strcspn(state->cache_line, " \t\r")
   ;  dim therank = 0
   ;  const char* theannot = ""
   ;  int id_in_query, flush_sequence

   ;  if (l < state->cache_length)
         state->cache_line[l]  =  '\0'
      ,  theannot =  state->cache_line+l+1       /* hackish, rework if necessary */

   ;  if (!l)
         miaow("empty identifier at line %u - ignoring sequence", state->linecount)
      ,  flush_sequence = 1
   ;  else 
         id_in_query    =  search_id(state->cache_line, thequery, &therank, state->linecount)
      ,  flush_sequence =  !id_in_query && universe      /* do not flush if found or subset */

   ;  if (flush_sequence)
      {  sylstatus status = SYL_SUCCESS
      ;  free(state->cache_line)
      ;  state->cache_line = NULL
      ;  state->n_skip++
                  /* fread_line will return NULL and set status to
                   * SYL_IGNORE when finding and ignoring sequence, it will
                   * return something and set status to SYL_SUCCESS when a
                   * new identifier is found, and return NULL and set status
                   * to SYL_DONE when EOF is found.
                  */
      ;  do
         {  int have_new_id = 0     /* currently not needed (fixme logic). */
         ;  dim linelen
         ;  state->cache_line
            =  fread_line(fp, &have_new_id, 1, &linelen, &(state->linecount), &status)
                  /* fixme: should fread_line, from cleanliness perspective
                   * not rather be tested on have_new_id here ?
                  */
         ;  if (status != SYL_IGNORE || state->cache_line)
            {  if (0) miaow("line %s %d", state->cache_line, (int) status)
            ;  break
         ;  }
                  /* With DONE cache_line can be != NULL (id without sequence at EOF).
                   * With DONE, sid == NULL is used as test in caller in case we
                   * need to ignore the last sequence - fixme.
                  */
      ;  }
         while (1)
      ;  return status ? status : SYL_IGNORE    /* so never return SYL_SUCCESS here */
   ;  }

      dest->sid   =  state->cache_line
   ;  state->cache_line = NULL
   ;  dest->annot =  theannot

   ;  if (id_in_query)
      dest->rank = therank
   ;  else
      dest->rank = state->rank_notfound++
   ;  return SYL_SUCCESS
;  }


sylstatus read_unit
(  INPUT       fp
,  unit   *    dest
,  query  *    thequery
,  int         universe       /* boolean */
,  struct read_state * state
)
   {  unsigned char* theseq =  NULL
   ;  dim seqlen        =  0
   ;  char dummy[]      =  "dummy"
   ;  void* mem_alloc   =  dummy       /* mem_alloc will test alloc/realloc success */
   ;  sylstatus status  =  SYL_SUCCESS
   ;  int have_new_id   =  0

   ;  dest->sid = NULL                 /* used as test for joint EOF and IGNORE */

   ;  if (!state->cache_line)
      return SYL_DONE

   ;  if (status)
      {  miaow("please test me, found status %d in read_unit", (int) status)
      ;  return status
   ;  }

      do
      {  status = read_id(fp, dest, thequery, universe, state)
   ;  }
      while (status == SYL_IGNORE)

   ;  if (status)          /* DONE or FAIL, IGNORE should be impossible  */
      return status        /* fixme: memclean ? */

   ;  if (!(mem_alloc = emptyseq()))
      return SYL_FAIL

   ;  theseq = mem_alloc

   ;  while (1)
      {  dim linelen = 0
      ;  char* line = fread_line(fp, &have_new_id, 0, &linelen, &(state->linecount), &status)

      ;  if (status == SYL_FAIL)    /* fixme: check/free line ? */
         break

                           /* noteme: we only set cache_line in the
                            * branch below, which garantuees SUCCESS on exit
                            * of this function.
                           */
      ;  if (!line || have_new_id)
         {  state->cache_line =  line
         ;  state->cache_length  =  linelen

         ;  dest->bases       =  theseq
         ;  dest->n_bases     =  seqlen
         ;  break
      ;  }

         {  if (!(mem_alloc = myrealloc(theseq, seqlen+linelen+1)))
            {  status = SYL_FAIL
            ;  break
         ;  }
            theseq = mem_alloc
         ;  memcpy(theseq+seqlen, line, linelen)
         ;  theseq[seqlen+linelen] = '\0'
         ;  seqlen += linelen
         ;  free(line)
         ;  line = NULL
      ;  }
      }
      return status
;  }


sylstatus read_fasta
(  collection* cl
,  INPUT       fp
,  query*      thequery
,  dim*        rank_max
,  unsigned    modes
,  int         length_co
,  int         max_cycle_length
,  int         gap_size_required
,  int         theK
)
   {  dim list_n_alloc  =  4        /* number of units alloced. cl->list_n tracks #used */
   ;  dim rm            =  0
   ;  dim n_rm          =  0
   ;  char dummy[]      =  "dummy"
   ;  void* mem_alloc   =  dummy    /* mem_alloc will test alloc/realloc success */
   ;  sylstatus status  =  SYL_SUCCESS
   ;  unit  receiver    =  { 0 }
   ;  int have_new_id   =  0
   ;  struct read_state state =  { NULL, 0, 0, 0, 0}
   ;  int warn_funny    =  modes & MODE_WARN_FUNNY
   ;  int universe      =  modes & MODE_UNIVERSE
   ;  struct mask_stats ms
   ;  int i

   ;  if ((fnrepeat_g = getenv("SYL_FN_REPEAT")))
      {  if (!strcmp(fnrepeat_g, "-"))
         fprepeat_g = stdout
      ;  else
         fprepeat_g = myfopen(fnrepeat_g, "w")
   ;  }

      state.rank_notfound = thequery->qel_n_file

   ;  for (i=0;i<16;i++)
         ms.n_nucs[i] = 0
      ,  ms.n_seqs[i] = 0
      ,  ms.n_pats[i] = 0

   ;  cl->list_n        =  0        /* number of sequences used */

   ;  state.cache_line
      =  fread_line(fp, &have_new_id, 0, &(state.cache_length), &(state.linecount), &status)

                                    /* v fixme: cache_line not freed in various paths below */
   ;  if (status == SYL_DONE)       /* empty FASTA file is dealt with at higher level */
      return SYL_SUCCESS
   ;  else if (status)
      return status

   ;  if (state.cache_line && !have_new_id)             /* clean-up: flush empty lines and encapsulate? */
      {  miaow("expect identifier at first line")
      ;  free(state.cache_line)
      ;  return SYL_FAIL
   ;  }

      if (!(cl->list = mymalloc(list_n_alloc * sizeof cl->list[0])))
      return SYL_FAIL

   ;  while (status != SYL_DONE)
      {  dim rmed = 0
      ;  status = read_unit(fp, &receiver, thequery, universe, &state)
      ;  if (status == SYL_FAIL || (status == SYL_DONE && !receiver.sid))
         break
                        /* ^ conceivably, we are DONE but want to ignore.
                         * this interface is somewhat ugly.
                        */

      ;  if (cl->list_n >= list_n_alloc)
         {  list_n_alloc *= 1.415
         ;  if (!(mem_alloc = myrealloc(cl->list, list_n_alloc * sizeof cl->list[0])))
            {  status = SYL_FAIL
            ;  break
         ;  }
            cl->list = mem_alloc
      ;  }
         memcpy(cl->list+cl->list_n, &receiver, sizeof cl->list[0])
      ;  rmed = validate_sequence
         (  cl->list+cl->list_n
         ,  warn_funny
         ,  length_co
         ,  max_cycle_length
         ,  gap_size_required
         ,  theK
         ,  modes
         ,  &ms
         )
      ;  if (rmed)
            rm += rmed
         ,  n_rm++
      ;  cl->list_n++
   ;  }

      if (status == SYL_FAIL)
      {  collection_release(cl)
      ;  if (state.cache_line)
            free(state.cache_line)
         ,  miaow("freeeeeing")
   ;  }

      else
      {  tell(VERBOSE_ONCE, "read %u sequences", (unsigned) cl->list_n)
      ;  if (thequery->qel_list)
         tell
         (  VERBOSE_ONCE
         ,  "searched/found/skipped %u/%u/%u ids"
         ,  (unsigned) thequery->qel_n_uniq
         ,  (unsigned) thequery->qel_n_found
         ,  (unsigned) state.n_skip
         )
      ;  if (max_cycle_length)
         {  tell
            (  VERBOSE_ONCE
            ,  "repeat-masked %lu bases in %lu sequences in %s mode"
            ,  (ulong) rm
            ,  (ulong) n_rm
            ,  gap_size_required ? "gappy" : "regular"
            )
         ;  tell
            (  VERBOSE_ONCE
            ,  " R   #seq    #frag      #nucl"
            )
         ;  for (i=1;i<=max_cycle_length;i++)
            tell
            (  VERBOSE_ONCE
            ,  "%2d %6d %8d %10d"
            ,  (int) i
            ,  (int) ms.n_seqs[i]
            ,  (int) ms.n_pats[i]
            ,  (int) ms.n_nucs[i]
            )
      ;  }

         *rank_max = state.rank_notfound
      ;  status = SYL_SUCCESS
   ;  }

      if (fprepeat_g)
      myfclose(&fprepeat_g)
   ;  return status
;  }


