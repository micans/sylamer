

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


int main(int argc, char* argv[]) {

   int help = argc == 2 && !strcmp(argv[1], "-h");
   FILE *fpin = myfopen(argc > 1 && !help ? argv[1] : "-", "r");
   char kmer[4096];
   unsigned linenum = 0;

   if (help) {
      fprintf(stderr, "Stream format:  LABEL  n1  n1+n2  k1  sample\n");
      exit(0);
   }

   while (!feof(fpin)) {
      unsigned n1, bg, k, sample;
      double hyperQ = 1.0, hyperP = 1.0, binomQ = 1.0, binomP = 1.0;
      linenum++;

      if (5 != fscanf(fpin, "%4095s%u%u%u%u\n", kmer, &n1, &bg, &k, &sample)) {
         if (0 && feof(fpin))
         break;
         die("could not parse line %u", linenum);
      }

      if (!bg)          die   ("zero bg (line %u)"    ,  linenum);
      if (n1 > bg)      die   ("n1 > bg (line %u)"    ,  linenum);
      if (sample > bg)  die   ("sample > bg (line %u)",  linenum);
      if (k > sample)   die   ("k > sample (line %u)" ,  linenum);
      if (k > n1)       die   ("k > n1 (line %u)"     ,  linenum);

      if (k)
      {  hyperQ = gsl_cdf_hypergeometric_Q(k-1, n1, bg-n1, sample);
         binomQ = gsl_cdf_binomial_Q(k-1, n1 * 1.0 / bg, sample);
      }

      if (k < sample)
      {  hyperP = gsl_cdf_hypergeometric_P(k, n1, bg-n1, sample);
         binomP = gsl_cdf_binomial_P(k, n1 * 1.0 / bg, sample);
      }

      fprintf
      (  stdout
      ,  "%s\t%s\t%g\t%g\t%g\t%g\t%u\t%u\t%u\t%u\n"
      ,  kmer
      ,  hyperP < hyperQ ? "ur" : "or"
      ,  hyperP
      ,  hyperQ
      ,  binomP
      ,  binomQ
      ,  (unsigned) n1
      ,  (unsigned) bg
      ,  (unsigned) k
      ,  (unsigned) sample
      );
   }
   return 0;
}




