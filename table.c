
#include <stdio.h>
#include <float.h>

#include "table.h"
#include "util.h"

static void heap_min_insert
(  float* heap
,  dim    capacity
,  float  val
,  dim*   n_inserted
)
   {  if (*n_inserted < capacity)
      {  dim d = *n_inserted

      ;  while (d != 0 && *(heap+(d-1)/2) < val)
         { *(heap+d) = *(heap+(d-1)/2)
         ;  d = (d-1)/2
      ;  }
         *(heap+d) = val
      ;  (*n_inserted)++
   ;  }
      else if (val < *heap)
      {  dim root = 0
      ;  dim d

      ;  while((d = 2*root+1) < capacity)
         {  if (*(heap+d) < *(heap+d+1))
            d++
         ;  if (val < *(heap+d))
            {  *(heap+root) = *(heap+d)
            ;  root = d
         ;  }
            else
            break
      ;  }
         *(heap+root) = val
   ;  }
   }


static void heap_max_insert
(  float* heap
,  dim    capacity
,  float  val
,  dim*   n_inserted
)
   {  if (*n_inserted < capacity)
      {  dim d =  *n_inserted

      ;  while (d != 0 && *(heap+(d-1)/2) > val)
         { *(heap+d) =  *(heap+(d-1)/2)
         ;  d = (d-1)/2
      ;  }
         *(heap+d) =  val
      ;  (*n_inserted)++
   ;  }
      else if (val > *heap)
      {  dim root  =  0
      ;  dim d

      ;  while((d = 2*root+1) < capacity)
         {  if (*(heap+d) > *(heap+d+1))
            d++
         ;  if (val > *(heap+d))
            {  *(heap+root) = *(heap+d)
            ;  root = d
         ;  }
            else
            break
      ;  }
         *(heap+root) = val
   ;  }
;  }


unsigned char* get_toptable
(  float* table
,  dim    n_rows
,  dim    n_cols
,  dim    n_cols_used
,  dim    top
)
   {  dim i, j
   ;  unsigned char* select = mycalloc(n_rows, sizeof select[0])

   ;  float* heap_min=  mycalloc(top+1, sizeof heap_min[0])    /* need 1 extra if top is even (fixme-why?) */
   ;  float* heap_max=  mycalloc(top+1, sizeof heap_max[0])    /* need 1 extra if top is even */
   ;  float* row_max =  mycalloc(n_rows, sizeof row_max[0])
   ;  float* row_min =  mycalloc(n_rows, sizeof row_min[0])
   ;  float  all_max =  -FLT_MAX
   ;  float  all_min =  FLT_MAX

   ;  dim n_min_inserted = 0
   ;  dim n_max_inserted = 0

   ;  if (!select || !heap_min || !heap_max || !row_max || !row_min)
      return NULL          /* we leak memory unless all callocs failed */

   ;  if (!(top & 1))      /* it's even .. */
         heap_min[top] = -FLT_MAX               /* need dummy elements */
      ,  heap_max[top] =  FLT_MAX               /* need dummy elements */

   ;  for (i=0;i<n_rows;i++)
      {  float this_max = -FLT_MAX
      ;  float this_min =  FLT_MAX

      ;  for (j=i*n_cols;j<i*n_cols+n_cols_used;j++)
         {  float f = table[j]
         ;  if (f > this_max) this_max = f
         ;  if (f < this_min) this_min = f
      ;  }

         heap_min_insert(heap_min, top, this_min, &n_min_inserted)
      ;  heap_max_insert(heap_max, top, this_max, &n_max_inserted)

      ;  row_max[i] = this_max
      ;  row_min[i] = this_min

      ;  if (all_max < this_max)
         all_max = this_max
      ;  if (all_min > this_min)
         all_min = this_min
   ;  }

      {  float thr_min = heap_min[0]
      ;  float thr_max = heap_max[0]
      ;  dim n_min_selected = 0
      ;  dim n_max_selected = 0
      ;  tell(VERBOSE_ONCE, "top-%d thresholds are %g and %g", (int) top, thr_max, thr_min)
      ;  tell(VERBOSE_ONCE, "extrema are %g and %g", all_max, all_min)
      ;  free(heap_min)
      ;  free(heap_max)
      ;  for (i=0;i<n_rows;i++)
         {  if (row_max[i] >= thr_max && n_max_selected < top)
               select[i] |= 2
            ,  n_max_selected++
         ;  if (row_min[i] <= thr_min && n_min_selected < top)
               select[i] |= 1
            ,  n_min_selected++
      ;  }
         free(row_max)
      ;  free(row_min)
   ;  }
#if 0                               /* don't ask. */
      frobotz(1.0, "table")
#endif                              /* (nonsense) */
   ;  return select
;  }



