
#ifndef  include_query
#define  include_query

#include "inttypes.h"
#include "sylio.h"

                        /* query element */
typedef struct
{  char* sid            /* sequence ID */
;  ofs   rank           /* occurrence */
;  int   found          /* found in fasta file? */
;
}  qel   ;


typedef struct
{  qel*  qel_list
;  dim   qel_n_uniq
;  dim   qel_n_file
;  dim   qel_n_found
;
}  query ;


int qel_cmp_by_id(const void* p1, const void* p2);
int qel_cmp_by_rank(const void* p1, const void* p2);
sylstatus read_query(INPUT fp, query* qry, int reverse);

void query_init(query* q);
void query_cat(query* q);
void query_release(query* q);
void query_uniq(query* q);


#endif

