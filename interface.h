
#ifndef  include_interface
#define  include_interface


#include "inttypes.h"
#include "sylio.h"
#include "unit.h"
#include "array.h"
#include "acc.h"
#include "words.h"


   /* As a general rule, nothing in this struct should be malloce'd.
    * anything on the heap should go in the data frame.
   */
struct parameter
{  unsigned long  modes
;  unsigned long  modes2
#define     MODE_WARN_FUNNY   1 <<  0        /* warn on bases other than A C G T U N X */
#define     MODE_R2ON         1 <<  1        /* erase 1-repeats and 2-repeats          */
#define     MODE_R2USERON     1 <<  2        /* do the above even if not default action*/
#define     MODE_LOGFOLD      1 <<  3        /* output log ratio scores                */
#define     MODE_REVERSE      1 <<  4        /* reverse the universe                   */
#define     MODE_2STRAND      1 <<  5        /* consider both strands, double sequences*/
#define     MODE_2STRAND2     1 <<  6        /* same, append reverse complement        */
#define     MODE_SHUFFLE      1 <<  7        /* shuffle sequences                      */
#define     MODE_UNIVERSE     1 <<  8        /* universe mode if on, subset otherwise  */
#define     MODE_COLNAMES     1 <<  9        /* output colnames                        */
#define     MODE_BYSTRETCH    1 << 10        /* count universe as non-overlapping kmers*/
#define     MODE_ULTERIOR     1 << 11        /* keep ulterior parts of UTR             */
#define     MODE_NONE         1 << 12        /* do not output anything (all modes? docme!) */
#define     MODE_FREAD        1 << 13        /* use fread buffering (block reads)      */
#define     MODE_PRINTID      1 << 14        /* print sequence IDs                     */
#define     MODE_STOP_READ    1 << 15        /* stop after reading data (timing test)  */
#define     MODE_STOP_BG      1 << 16        /* stop after computing background (idem) */
#define     MODE_STOP_WORDS   1 << 17        /* stop after computing and writing words */
#define     MODE_OVER         1 << 18        /* output overrepresentation              */
#define     MODE_UNDER        1 << 19        /* output underrepresentation             */
#define     MODE_BINOMIAL     1 << 20        /* compute binomial P-values (long-listing)*/
#define     MODE_NOANCHOR     1 << 21        /* do not anchor markov correction        */
#define     MODE_FREQUENCIES  1 << 22        /* read frequencies along wordlist        */
#define     MODE_AND          1 << 23        /* all conditions be met                  */
#define     MODE_SUM          1 << 24        /* all conditions be met                  */
#define     MODE_BAG          1 << 25        /* baggify all words                      */
#define     MODE_FSTAT        1 << 26        /* summarize fasta file                   */
#define     MODE_PERTURB_TIES 1 << 27        /*                                        */
#define     MODE_WSTAT2       1 << 28        /*                                        */
#define     MODE_GAPPY        1 << 29        /* are we in gapped (== dyad) mode?       */
#define     MODE_MASKALL      1 << 30        /* mask repeats fully                     */
#define     MODE_SORTBYLENGTH 1 << 31        /* sort by length!                        */

#define     Mode_BINBYLENGTH  1 <<  0        /* bin by length .                        */  
#define     Mode_UTRATTR      1 <<  1        /* add utr attributes to output table     */
#define     Mode_STEPPING     1 <<  2        /* do not nest bins, do consecutive bins  */
;  unsigned long dyad_modes
#define     DYAD_REPEAT       1 << 0
#define     DYAD_INVERT       1 << 1
#define     DYAD_UNIFORM      1 << 2         /* use expected on uniform dyad component distribution */
;  unsigned long main_modes
#define     MAIN_EXTREME      1 <<  0        /* output max/min                         */
#define     MAIN_TABLE        1 <<  1        /* provide table output                   */
#define     MAIN_FILTER       1 <<  2        /* one of the main sylamer modes          */
#define     MAIN_COUNT        1 <<  3        /* one of the main sylamer modes          */
#define     MAIN_WSTAT        1 <<  4        /* produce word frequencies               */
#define     MAIN_MWW          1 <<  5        /* Mann-Whitney-Wilcoxon test             */
;  double      cutoff
;  unsigned    per_sequence         /* count #utrs with >= per_sequence_p hit (only if set) */
;  unsigned    cap                  /* allow now no more than cap word-utr hits (only if set) */
;  unsigned    lastbin
;  const char* format
;  int         unit_at_least
;  int         unit_at_most
;  unsigned    shift_left           /* number of leading things to ignore  */
;  unsigned    shift_right          /* number of trailing things to ignore */

;  unsigned    binsize
;  unsigned    splitsize
;  unsigned    split_remainder

#if 0
;  unsigned    accelerate_num       /* increment bin sizes this many times */
;  unsigned    accelerate_offset    /* initial segment size with constant bin size */
;  unsigned    accelerate_remainder
;  unsigned    accelerate_rampsize
#endif

;  unsigned    n_universe_fake      /* fake universe size */
;  unsigned    utimes_fake          /* fake universe size as a multiple of set size */
;  unsigned    unit_size_fake       /* fake unit size */
;  int         length_co            /* cut off sequences longer than this */
;  int         repeat_mask          /* mask repeats of this wave-length (e.g. dimer, monomer) */

;  dim         toptable
;  dim         n_trial

;  unsigned    K_g                  /* The k in kmer. It rules them all. _g just to make it look bigger.  */
;  dim         M_g                  /* Markov word length (order + 1)      */

;  dim         gap_min
;  dim         gap_max
;  dim         gap_gap

;  FILE        *iffplog 
;  FILE        *fpresult

;  INPUT       fpfasta
;  INPUT       fpquery
;  INPUT       fpwords

;  const char* fnwords
;  const char* fnresult
;  const char* fnfound
;  const char* fnwquery
;  const char* fnmissing
;  const char* fnexp

;  const char* patlist
;  const char* xpatlist
;  double entropy

;  const char* test
;  dim hash_commit_fac_g
;  dim hash_commit_min_g
; 
}  ;

extern double out_of_bounds_g;

struct bininfo
{  dim      offset
;  dim      size
;  double   utr_length
;  double   gccontent
;  double   n_masked
;
}  ;


struct dataframe
{  collection* thecollection
;  dim collection_n_found           /* number of things found in query                       */
                           /* the data layout and handling for thecollection is not
                            * in a very tight and rigid and easily readable state yet.
                            * This pertains to how it is sorted, how the
                            * ranks are computed, and where the rank threshold
                            * is stored (that is, in collection_n_found)
                           */
;  dim rank_max                     /* one less than maximum assigned rank in read_fasta     */
                                    /* -> it's presence in this struct is somewhat doubtful  */
;  query*  thequery

;  harray* hash_template_g
;  unsigned hmask_g                 /* specifies integer hash array size   */

;  float* eoracle_g                 /* storage for expected frequencies */
;  wlist wordlist

;  float* extrema
;  float* table                     /* for log transformed P-values */
;  float* utrtable                  /* for MWW mode, need value for all UTRs and all words */
;  struct bininfo*  bins_info       /* bin summary statistics */

;  dim*  tablect                    /* for counts */
;  dim   bg_use_count               /* the number to use for background counts; depends on modality */
;  dim   n_analyse                  /* number of sequences to analyse */
;  dim   n_universe                 /* number of sequences in universe */
;  dim   n_bin                      /* number of bins as precomputed (could be overridden by lastbin) */

;  acc   bg
;  acc   bg_gap_cache
;  acc   bg_xy
;  acc   bg_xyz
;  acc   window
;  acc   window_gap_cache
;  acc   window_xy
;  acc   window_xyz
;  mkv   mkv1
;  mkv   mkv2

;  struct  parameter* par
;
}  ;


         /* nothing is alloc'ed, everything is set to 0 or NULL */
void dataframe_null
(  struct dataframe* dt
,  struct parameter* par
,  collection *cl
,  query *qry
)  ;

int dataframe_release(struct dataframe* dt);


         /* this struct encodes both user settings as
          * well as session settings related to binning.
         */
void parameter_init(struct parameter* par);



unsigned parse_args
(  struct parameter* par
,  int argc
,  char* argv[]
)  ;



#endif

