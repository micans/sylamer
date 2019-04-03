
#include "interval.h"
#include "util.h"

unsigned interval_prepare_bywidth
(  collection* thecl
,  dim n_analyse
,  unsigned binsize
)
   {  dim i, n_bin=0
               /* n_analyse 10,050     n_analyse 120
                * binsize 200          binsize   200
                * rem 50               rem       120
               */
   ;  dim rem = n_analyse > binsize ? n_analyse % binsize : 0
   ;  for (i=0;i<n_analyse;i++)
      {  dim delta = i > n_analyse / 2 ? rem : 0
      ;  if ((i-delta) % binsize == 0)
            thecl->list[i].firstinbin = 1
         ,  n_bin++
   ;  }
      return n_bin
;  }


unsigned interval_prepare_bycount
(  collection* thecl
,  dim n_analyse
,  unsigned splitsize
,  int binbylength
)
   {  dim i = 0, totalsize = 0, n_window, size_of_window
   ;  double average_chunk_size
   ;  for (i=0;i<n_analyse;i++)
      {  unit* ut = thecl->list+i
      ;  totalsize += binbylength ? ut->n_bases : 1
   ;  }
   ;  if (!totalsize)
      return 1
   ;  average_chunk_size = totalsize * 1.0 / splitsize
 
   ;  size_of_window = binbylength ? thecl->list[0].n_bases : 1
   ;  thecl->list[0].firstinbin = 1
   ;  n_window = 1

   ;  for (i=1;i<n_analyse;i++)
      {  int lookahead = binbylength && i+1 < n_analyse
      ;  double diffthis, diffnext
      ;  size_of_window += binbylength ? thecl->list[i].n_bases : 1
      ;  diffthis = n_window * average_chunk_size - size_of_window
      ;  diffnext =     lookahead
                     ?     ( size_of_window + 1.0 * thecl->list[i+1].n_bases )
                        -  ( n_window * average_chunk_size )
                     :  0.0
      ;  if (diffthis < 0 || (diffnext > diffthis))
         {  thecl->list[i].firstinbin = 1
;if(0)fprintf(stderr, "bin start %d - %d diff %g %g\n", (int) i, (int) size_of_window,  diffthis, diffnext)
         ;  n_window++
      ;  }
      }
      return n_window
;  }

unsigned interval_prepare
(  collection* thecl
,  dim n_analyse
,  unsigned splitsize
,  unsigned binsize
,  int binbylength
)
   {  if (splitsize)
      return interval_prepare_bycount(thecl, n_analyse, splitsize, binbylength)
   ;  else
      return interval_prepare_bywidth(thecl, n_analyse, binsize)
;  }


void interval_step
(  collection* thecl
,  dim n_analyse
,  spacer* sp
)
   {  unit* ut = thecl->list + sp->i_hi + 1

   ;  if (sp->i_hi)
      sp->i_lo = sp->i_hi

   ;  while (ut < thecl->list+n_analyse && !ut->firstinbin)
      ut++
   ;  sp->i_hi = ut - thecl->list
;  }


