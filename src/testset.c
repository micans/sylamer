

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include "gsl/gsl_cdf.h"


void die
(  const char* fmt
,  ...
)
  {  va_list  args
  ;  va_start(args, fmt)
  ;  vfprintf(stderr, fmt, args)
  ;  fputc('\n', stderr)
  ;  va_end(args)
  ;  exit(1)
; }


FILE* myfopen (const char* fname, const char* mode) {
   FILE* fp = NULL;
   if (!strcmp(fname, "-")) {
      if (!strcmp(mode, "w")) {
         fp = stdout;
      }
      else if (!strcmp(mode, "r")) {
         fp = stdin;
      }
      else {
         die("fire programmer: unrecognized open mode <%s>", mode);
      }
   }
   else {
      fp = fopen(fname, mode);
      if (!fp) die("cannot open file <%s> in mode <%s>", fname, mode);
   }
   return fp;
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



int main(int argc, char* argv[]) {

#define N_MAX 30000
   rank_unit *rus = malloc(N_MAX * sizeof rus[0]);
   char buf[1000];
   ulong n_rus = 0, i = 0;
   ulong stretch;

   while (fgets(buf, sizeof buf, stdin)) {
      if (n_rus > N_MAX) {
         die("buffer size exceeded");
      }
      rus[n_rus].index = n_rus;
      rus[n_rus].value = atof(buf); 
      rus[n_rus].ord = -1.0;
      n_rus++;
   }

   qsort(rus, n_rus, sizeof rus[0], rank_unit_cmp_value);

   rus[n_rus].value = rus[n_rus-1].value + 1000;    /* sentinel out-of-bounds value */
   stretch = 1;
   for (i=0;i<n_rus;i++) {
      if (rus[i].value != rus[i+1].value) {         /* input current stretch */
         unsigned j;
         for (j=0;j<stretch;j++) {
            rus[i-j].ord = 1 + (i+1 + i-stretch) / 2.0;
         }
         stretch = 1;
      }
      else
      stretch++;
   }
   qsort(rus, n_rus, sizeof rus[0], rank_unit_cmp_index);

   {  long n1 = 0, n2 = n_rus;
      double R1 = 0.0, U = 0.0, sigma, mean, p;
      for (i=0;i<n_rus-1;i++) {
         int shiftleft;
         n1++; n2--;
         mean = n1 * (n2) * 0.5;
         sigma = sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0);
         R1 += rus[i].ord;
         U = R1 - (n1 * (n1 + 1.0)) * 0.5;
         shiftleft = U < mean;
         p =   U < mean
            ?  gsl_cdf_ugaussian_P(( U - mean ) / sigma)
            :  gsl_cdf_ugaussian_Q(( U - mean ) / sigma);
         fprintf(stdout, "%d\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\n", (int) i, rus[i].value, rus[i].ord, U, R1, p, (shiftleft ? -1 : 1 ) *  log(p) / log(10));
      }
   }
   return 0;
}


/*                 n1 * (n1 + 1)
 *    U  =   R1 -  -------------
 *                      2

 *    U - m          n1 n2            /   n1 n2 ( n1 + n2 + 1)
 *    -----    m  =  -----   s  = \  /    --------------------
 *      s              2           \/              12

*/



