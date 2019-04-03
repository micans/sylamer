
#ifndef  include_array
#define  include_array

#include "inttypes.h"
#include "words.h"


typedef union
{  unsigned num
;  float    val
;
}  rcell    ;


typedef struct
{  unsigned hid      /* the numerical representation of the kmer */
;  rcell    pl       /* payload */
;
}  hcell    ;

                     /* hash array. See implementation notes. Note: words is all
                      * the words represented in the hash. Ownership lies *not*
                      * with the hash (each hash has the same shallow copy).
                      * words determines iteration order.
                     */
typedef struct
{  hcell*   cells
;  unsigned hmask
;  const unsigned* words
;  dim      n_words
;
}  harray   ;


                     /* Search hash for a word. sets offset either to
                      * the slot where the word was found, or to
                      * an empty slot (indicated by hid == -1u)
                     */
#define H_ITERATE(thehash, thehid)              \
         hcell* cells =  thehash->cells         \
      ;  unsigned mask=  thehash->hmask         \
      ;  unsigned hhh =  hashwell(thehid,mask)  \
      ;  dim jj = 0                             \
      ;  while(cells[hhh].hid != -1u && cells[hhh].hid != thehid) \
         {  hhh += HASH_JUMP                    \
         ;  hhh &= mask                         \
         ;  if (jj++ > mask) return arrr("hash error") \
      ;  }

#define HSET(thehash, thehid, thenum)           \
      {  H_ITERATE(thehash, thehid)             \
         if (cells[hhh].hid == -1u)             \
         cells[hhh].hid = thehid                \
      ;  cells[hhh].pl.num = thenum             \
   ;  }

#define HADD_IF_PRESENT(thehash, thehid, inc)   \
      {  H_ITERATE(thehash, thehid)             \
         if (cells[hhh].hid == thehid)          \
         cells[hhh].pl.num += inc               \
   ;  }
#define HINC_IF_PRESENT(thehash, thehid)     HADD_IF_PRESENT(thehash, thehid, 1)

         /* This macro retrieves a value from a hash cell
          * and sets the hash cell value to thenum.
         */
#define HSET_IF_PRESENT(lvalue, thehash, thehid, thenum)   \
      {  H_ITERATE(thehash, thehid)             \
      ;  lvalue = 0                             \
      ;  if (cells[hhh].hid == thehid)          \
         {  lvalue = cells[hhh].pl.num          \
         ;  cells[hhh].pl.num = thenum          \
      ;  }                                      \
      }
#define HZERO_IF_PRESENT(lvalue, thehash, thehid)   HSET_IF_PRESENT(lvalue, thehash, thehid, 0)

#define HOFFSET_IF_PRESENT(lvalue, thehash, thehid)        \
      {  H_ITERATE(thehash, thehid)             \
         if (cells[hhh].hid == thehid)          \
         lvalue = hhh                           \
   ;  }

#define HISPRESENT(lvalue, thehash, thehid)     \
      {  H_ITERATE(thehash, thehid)             \
      ;  lvalue = cells[hhh].hid == thehid      \
   ;  }


                     /* regular array */
typedef struct
{  rcell*   cells
;  dim      n_cells
;
}  rarray   ;


         /* Expected occurrence frequencies.
         */
typedef struct
{  rarray  *rval
;  harray  *hval
;
}  mkv      ;


rarray* rarray_new(unsigned k);
void rarray_reset (rarray* ra);
void rarray_free(rarray* ra);

harray* harray_new(wlist* wl, unsigned mask, unsigned K);
harray* harray_clone (harray* ha);
void harray_reset (harray* ha);
void harray_free(harray* ha);

ofs array_step
(  ofs         offset
,  harray*     ha
,  rarray*     ra
,  unsigned*   hidp
,  unsigned*   hhhp        /* the pivotal hash index, if applicable */
,  dim*        countp
,  float**     valpp
,  sylstatus   *mystatus
)  ;

void mkv_null
(  mkv* markov
)  ;

sylstatus mkv_init
(  mkv* markov
,  unsigned k
,  harray* ha
)  ;

void mkv_free
(  mkv* markov
)  ;

ofs mkv_iterate
(  ofs         offset      /* for either markov->hval or markov->rval */
,  mkv*        markov
,  unsigned*   hidp        /* store the current word                  */
,  float**     valpp       /* make available for writing              */
,  sylstatus*  mystatus
)  ;


#endif

