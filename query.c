
#include "query.h"
#include "sylio.h"
#include "util.h"

#include <string.h>
#include <stdio.h>


void query_init
(  query* thequery
)
   {  thequery->qel_list     =  NULL
   ;  thequery->qel_n_uniq   =  0
   ;  thequery->qel_n_file   =  0
   ;  thequery->qel_n_found  =  0
;  }


void query_release
(  query* thequery
)
   {  dim i
   ;  if (thequery)
      {  for (i=0;i<thequery->qel_n_file;i++)
         free(thequery->qel_list[i].sid)
      ;  free(thequery->qel_list)
   ;  }
      thequery->qel_n_file = 0
   ;  thequery->qel_list = NULL
;  }


int qel_cmp_by_id(const void* p1, const void* p2)
   {  const qel *qel1 = p1, *qel2 = p2
   ;  return strcmp(qel1->sid, qel2->sid)
;  }


int qel_cmp_by_rank(const void* p1, const void* p2)
   {  const qel *qel1 = p1, *qel2 = p2
   ;  return qel1->rank < qel2->rank ? -1 : qel1->rank == qel2->rank ? 0: 1
;  }


sylstatus read_query(INPUT fp, query* thequery, int reverse)
   {  char* line        =  NULL
   ;  dim linelen       =  0
   ;  dim linecount     =  0
   ;  dim qry_n_alloc   =  1000
   ;  dim qry_n         =  0
   ;  qel* qry          =  mymalloc(qry_n_alloc * sizeof qry[0])
   ;  qel* qry_realloc  =  qry
   ;  sylstatus status  =  SYL_SUCCESS

   ;  if (!qry)
      return arrr("could not alloc query")

   ;  while ((line = fread_line(fp, NULL, 0, &linelen, &linecount, &status)))
      {  size_t l = strcspn(line, " \t\r")

      ;  if (status && status != SYL_DONE)
         break

;if(0)fprintf(stderr, "status %d line [%s]\n", (int) status, line)

                  /*  fread_line should have skipped the empty line */
      ;  if (strlen(line) == 0)
         {  miaow("empty line while reading query -- (odd, my servant should skip those)")
         ;  free(line)
         ;  continue
      ;  }
         line[l] = '\0'
      ;  if (qry_n >= qry_n_alloc)
         {  qry_n_alloc *= 1.415       /* my very own favourite alloc constant */
         ;  qry_realloc = myrealloc(qry, qry_n_alloc * sizeof qry[0])
         ;  if (!qry_realloc)
            break
         ;  qry = qry_realloc
      ;  }
         qry[qry_n].sid = line
      ;  qry[qry_n].rank = qry_n
      ;  qry[qry_n].found = 0
      ;  qry_n++
      ;  line = NULL
   ;  }

      thequery->qel_n_file = qry_n
   ;  thequery->qel_list = qry

   ;  if (!qry_realloc || status != SYL_DONE)
      {  if (line)
         free(line)
      ;  query_release(thequery)
      ;  status = SYL_FAIL
   ;  }
      else
      {  tell(VERBOSE_ONCE, "read %u query ids", (unsigned) qry_n)

      ;  if (reverse)
         {  dim i
         ;  for (i=0;i<qry_n;i++)
            qry[qry_n-i-1].rank = i
      ;  }
      }
      return status
;  }


void query_cat(query * thequery)
   {  dim i
   ;  for (i=0;i<thequery->qel_n_file;i++)
      {  qel* cur = thequery->qel_list+i
      ;  fprintf(stdout, "%s %d %d\n", cur->sid, (int) cur->rank, (int) cur->found)
   ;  }
   }


void query_uniq(query * thequery)
   {  qel* cur = thequery->qel_list+0, *max = cur+thequery->qel_n_file, *piv, tmp
   ;  thequery->qel_n_uniq = 0

               /* we need below number > 0 for qsort AND for deduplicating */
   ;  if (thequery->qel_n_file == 0)
      return

   ;  qsort
      (  thequery->qel_list
      ,  thequery->qel_n_file
      ,  sizeof thequery->qel_list[0]
      ,  qel_cmp_by_id
      )

   ;  for (piv=thequery->qel_list+1;piv<max;piv++)
      {  int diff = strcmp(cur->sid, piv->sid)

      ;  if (diff)               /* pivot is a new element */
         {  cur++
         ;  if (cur != piv)
            {  tmp = cur[0]
            ;  cur[0] = piv[0]   /* move pivot forwards                       */
            ;  piv[0] = tmp      /* move some duplicate to pivot's position   */
         ;  }
         }
         else           /* pivot same as latest kept element */
         {  if (cur->rank > piv->rank) 
            cur->rank = piv->rank
         ;  tell
            (  VERBOSE_BIN
            ,  "query element %s now rank %d"
            ,  cur[0].sid
            ,  (int) cur->rank
            )
      ;  }
      }
      thequery->qel_n_uniq = cur-thequery->qel_list+1
;  }



