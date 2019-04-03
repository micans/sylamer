
#ifndef  include_unit
#define  include_unit

#include "inttypes.h"
#include "query.h"


                  /* This encodes a stretch of bases, which will
                   * usually be sequence of some kind, for example a UTR. But
                   * we avoid 'seq' throughout this source when naming
                   * variables.
                  */
typedef struct
{  char*       sid
;  const char* annot          /* warning: never free this (it's piggy-backing on sid, for now) */
;  unsigned char* bases
;  dim         n_bases
;  dim         n_masked
;  double      gccontent      /* fraction */
;  dim         rank           /* determined by occurrence in query */
;  unsigned char firstinbin   /* this unit defines start of bin */
;
}  unit  ;

                  /* unimaginatively named. A collection of units/sequences. */
typedef struct
{  unit* list
;  dim list_n
;
}  collection ;

int unit_cmp_by_rank(const void* p1, const void* p2);
int unit_cmp_by_revrank(const void* p1, const void* p2);
int unit_cmp_by_length(const void* p1, const void* p2);
int unit_cmp_by_revlength(const void* p1, const void* p2);

void collection_release (collection* cl);
sylstatus collection_revcompl(collection* cl);                 /* interweaves new sequences, doubles sequence count */
sylstatus collection_revcompl2(collection* cl, dim gap_max);   /* appends to sequences */
void collection_shuffle(collection* cl, dim n_universe);

sylstatus read_fasta
(  collection* cl
,  INPUT       fp
,  query*      thequery
,  dim*        rank_max
,  unsigned    modes
,  int         length_co
,  int         repeat_mask
,  int         gap_size_required
,  int         K
)  ;

   /* must be supplied with address of valid collection */
void collection_init
(  collection* cl
)  ;

#endif

