

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "interface.h"
#include "util.h"
#include "version.h"


double out_of_bounds_g = 100.0;


static unsigned syldone
(  unsigned status         /* SYL_DONE or SYL_FAIL */
,  const char* fmt
,  ...
)
   {  if (fmt)
      {  va_list  args
      ;  va_start(args, fmt)
      ;  vfprintf(stderr, fmt, args)
      ;  fputc('\n', stderr)
      ;  va_end(args)
   ;  }
      return status
;  }


         /* nothing is alloc'ed, everything is set to 0 or NULL */
void dataframe_null
(  struct dataframe* dt
,  struct parameter* par
,  collection *cl
,  query *qry
)
   {  wlist wl             =  { NULL, 0, 0 }
   ;  dt->thecollection    =   cl
   ;  dt->collection_n_found = 0
   ;  dt->rank_max         =   0

   ;  dt->thequery         =  qry

   ;  dt->hash_template_g  =  NULL
   ;  dt->hmask_g          =  0

   ;  dt->eoracle_g        =  NULL
   ;  dt->wordlist         =  wl

   ;  dt->extrema          =  NULL
   ;  dt->table            =  NULL
   ;  dt->utrtable         =  NULL
   ;  dt->bins_info        =  NULL
   ;  dt->tablect          =  NULL

   ;  dt->bg_use_count     =  0
   ;  dt->n_analyse        =  0
   ;  dt->n_universe       =  0
   ;  dt->n_bin            =  0

   ;  dt->par              =  par

   ;  acc_null(&(dt->bg))
   ;  acc_null(&(dt->bg_gap_cache))
   ;  acc_null(&(dt->bg_xy))
   ;  acc_null(&(dt->bg_xyz))
   ;  acc_null(&(dt->window))
   ;  acc_null(&(dt->window_gap_cache))
   ;  acc_null(&(dt->window_xy))
   ;  acc_null(&(dt->window_xyz))

   ;  mkv_null(&(dt->mkv1))
   ;  mkv_null(&(dt->mkv2))
;  }


int dataframe_release(struct dataframe* dt)
   {  if (dt->hash_template_g)
      harray_free(dt->hash_template_g)

   ;  acc_free(&(dt->bg))
   ;  acc_free(&(dt->window))
   ;  acc_free(&(dt->bg_gap_cache))
   ;  acc_free(&(dt->window_gap_cache))

   ;  if (!dt->par)
      miaow("dataframe release", "I seem to be called one time too often")
   ;  else
      {  if (dt->par->M_g)
         {  acc_free(&(dt->bg_xy))
         ;  acc_free(&(dt->bg_xyz))
         ;  acc_free(&(dt->window_xy))
         ;  acc_free(&(dt->window_xyz))
         ;  mkv_free(&(dt->mkv1))
         ;  mkv_free(&(dt->mkv2))
      ;  }

         if (dt->par->fpfasta) myfzclose(&(dt->par->fpfasta))
      ;  if (dt->par->fpquery) myfzclose(&(dt->par->fpquery))
      ;  if (dt->par->fpwords) myfzclose(&(dt->par->fpwords))
      ;  if (dt->par->fpresult) myfclose(&(dt->par->fpresult))
   ;  }

      query_release(dt->thequery)
   ;  collection_release(dt->thecollection)

   ;  wlist_free(&dt->wordlist)

   ;  if (dt->extrema)   free(dt->extrema)
   ;  if (dt->table)     free(dt->table)
   ;  if (dt->utrtable)  free(dt->utrtable)
   ;  if (dt->tablect)   free(dt->tablect)
   ;  if (dt->bins_info) free(dt->bins_info)
   ;  if (dt->eoracle_g) free(dt->eoracle_g)

   ;  dataframe_null(dt, NULL, NULL, NULL)
   ;  return 0
;  }


                           /* this struct encodes both user settings as
                            * well as session settings related to binning.
                           */
void parameter_init
(  struct parameter* par
)
   {  par->modes           =  MODE_R2ON | MODE_FREAD | MODE_BINOMIAL | MODE_WARN_FUNNY
   ;  par->modes2          =  Mode_BINBYLENGTH
   ;  par->main_modes      =  MAIN_TABLE
   ;  par->dyad_modes      =  0
   ;  par->cutoff          =  0.0
   ;  par->per_sequence    =  0
   ;  par->cap             =  0
   ;  par->lastbin         =  0
   ;  par->format          =  NULL
   ;  par->unit_at_least   =  -1
   ;  par->unit_at_most    =  -1
   ;  par->shift_left      =  0
   ;  par->shift_right     =  0

   ;  par->binsize         =  0
   ;  par->splitsize       =  0
   ;  par->split_remainder =  0

   ;  par->n_universe_fake =  0
   ;  par->utimes_fake     =  0
   ;  par->unit_size_fake  =  0
   ;  par->length_co       =  0
   ;  par->repeat_mask     =  0

   ;  par->toptable        =  0
   ;  par->n_trial         =  1

   ;  par->K_g             =  0 
   ;  par->M_g             =  0 

   ;  par->gap_min         =  1
   ;  par->gap_max         =  0
   ;  par->gap_gap         =  0

   ;  par->iffplog         =  NULL
   ;  par->fpresult        =  NULL

   ;  par->fpfasta         =  NULL
   ;  par->fpquery         =  NULL
   ;  par->fpwords         =  NULL

   ;  par->fnwords         =  NULL
   ;  par->fnresult        =  "-"
   ;  par->fnfound         =  NULL
   ;  par->fnwquery        =  NULL
   ;  par->fnmissing       =  NULL
   ;  par->fnexp           =  NULL

   ;  par->patlist         =  NULL     /* add these words */
   ;  par->xpatlist        =  NULL     /* and exclude those */
   ;  par->entropy         =  0.0

   ;  par->test            =  NULL
   ;  par->hash_commit_fac_g  =  4
   ;  par->hash_commit_min_g  =  2048
;  }



/* Lo and behold, a simple CPP command-line parsing framework
 * It employs these variables defined in sylmain:
 *    noargallowed, argv_i, argv, argc, dt.
 * It returns via mainfail() if an error is detected.
*/

#define thearg()  (  noargallowed      \
                  ?  miaow             \
                     (  "programmer brain damage, do not use option %s"    \
                     ,  argv[argv_i]   \
                     )                 \
                     ,  "argument-not-available" : argv[argv_i]     \
                  )
#define arg_switch()          while (++argv_i < argc) { if (0) {}
#define incopt(f)             if (++argv_i >= argc) return syldone(SYL_FAIL, "option %s requires argument", f)
#define optarg(f)             else if (!strcmp(thearg(), f)) { incopt(f);
#define uniarg(f)             else if (!strcmp(thearg(), f)) { noargallowed = 1;
#define uniargs(f,g)          else if (!strcmp(thearg(), f) || !strcmp(thearg(), g)) { noargallowed = 1;
#define failarg()             else { return syldone(SYL_FAIL, "unrecognized option <%s>", thearg()); }
#define endarg()              noargallowed = 0; }
#define arg_done()            }

#define set(mode)             (par->modes |= (mode))
#define unset(mode)           (par->modes |= (mode), par->modes ^= (mode))
#define set2(mode)            (par->modes2 |= (mode))
#define unset2(mode)          (par->modes2 |= (mode), par->modes2 ^= (mode))
#define main_set(mode)        (par->main_modes |= (mode))
#define main_unset(mode)      (par->main_modes |= (mode), par->main_modes ^= (mode))


static  const char* contributors[] =
{  "Anton Enright: design and use suggestions, suggestion to use C"
,  "Jacques van Helden: author of oligo-analysis, a related program"
,  "Weldon Whitener: bootstrap"
,  "Sergei Manakov: suggestions, testing"
,  "Harpreet Saini: testing"
,  "Russell Grocock: base feature set"
,  "The GSL authors/maintainers: high quality numerical computing"
};


unsigned parse_args
(  struct parameter* par
,  int argc
,  char* argv[]
)
   {  int argv_i = 0
   ;  int noargallowed = 0
   ;  int i

   ;  char* db = getenv("SYL_DB1")
   
   ;  syl_db1_G = db ? atof(db) : 0.0

   ;  arg_switch()

      uniarg("--no-anchor")      set(MODE_NOANCHOR);           endarg()

         uniarg("--stop-bg")     set(MODE_STOP_BG);            endarg()
         uniarg("--stop-read")   set(MODE_STOP_READ);          endarg()

         uniarg("--stop-words")
            set(MODE_STOP_WORDS);
            par->fpfasta = myfopen("/dev/null", "r");
         endarg()

      uniarg("--none")           set(MODE_NONE);               endarg()
      uniarg("--format-fasta")   par->format = ">%I%n%F%n";    endarg()
      uniarg("--format-id")      par->format = "%I%n";         endarg()
      optarg("-format")          par->format = thearg();       endarg()
      uniarg("--word-count")     main_set(MAIN_COUNT);         endarg()
      uniarg("--print-id")       set(MODE_PRINTID);            endarg()
      uniarg("--stat")           set(MODE_FSTAT);              endarg()
      uniarg("--word-stat")      main_set(MAIN_WSTAT);         endarg()
      uniarg("--word-utr-stat")  main_set(MAIN_MWW); set(MODE_WSTAT2); endarg()     /* WSTAT2 done in same block as MWW */
      uniarg("--mww")            main_set(MAIN_MWW);           endarg()
      uniarg("--under")          set(MODE_UNDER);              endarg()
      uniarg("--over")           set(MODE_OVER);               endarg()
      uniarg("--filter-and")     set(MODE_AND);                endarg()
      uniarg("--filter-sum")     set(MODE_SUM);                endarg()
      uniarg("--rm-maskall")     set(MODE_MASKALL);            endarg()

      uniarg("--alien-ok")       unset(MODE_WARN_FUNNY);       endarg()
      uniarg("--r2-off")         unset(MODE_R2ON);             endarg()
      uniarg("--r2-on")          set(MODE_R2USERON);           endarg()
      uniarg("--long-listing")   main_unset(MAIN_TABLE);       endarg()
      uniarg("--no-binomial")    unset(MODE_BINOMIAL);         endarg()

      uniarg("--ck")             set(MODE_BYSTRETCH);          endarg()
      uniarg("--utr-v")          set(MODE_ULTERIOR);          endarg()
      uniarg("--logfold")        set(MODE_LOGFOLD);            endarg()
      uniarg("--reverse")        set(MODE_REVERSE);            endarg()
      uniarg("--sort-by-length") set(MODE_SORTBYLENGTH);       endarg()
      uniarg("--bin-by-seqs")    unset2(Mode_BINBYLENGTH);     endarg()
      uniarg("--two-strands")    set(MODE_2STRAND);            endarg()
      uniarg("--twostrands")     set(MODE_2STRAND2);           endarg()
      uniarg("--shuffle")        set(MODE_SHUFFLE);            endarg()
      uniarg("--perturb-ties")   set(MODE_PERTURB_TIES);       endarg()
      uniarg("--utr-attributes") set2(Mode_UTRATTR);           endarg()

      uniarg("--dyad-repeat")    par->dyad_modes |= DYAD_REPEAT;  endarg()
      uniarg("--dyad-invert")    par->dyad_modes |= DYAD_INVERT;  endarg()
      uniarg("--dyad-paired")    par->dyad_modes |= DYAD_REPEAT | DYAD_INVERT;  endarg()
      uniarg("--monad-uniform")  par->dyad_modes |= DYAD_UNIFORM;  endarg()
      uniarg("--bag")            set(MODE_BAG);                endarg()

      uniarg("--version")
fprintf
   (  stdout,
      "Sylamer version %s, "
      "Copyright 2008 Wellcome Trust Sanger Institute\n",
      SYLAMER_VERSION
   );
puts("You may redistribute copies of sylamer under the terms of the GNU");
puts("General Public License, version 3 or later.");
puts("Authors:\n   Stijn van Dongen and Cei Abreu-Goodger");
puts("Contributors:");
for(i=0;i<sizeof contributors/sizeof contributors[0];i++) {
fputs("   ", stdout); puts(contributors[i]); }
return SYL_DONE;
      endarg()

      uniarg("--format-syntax")
puts("This documents the syntax for the -format option");
puts("%I       identifier");
puts("%S       sequence (on a single line - useful for streaming)");
puts("%F       formatted sequence (60 bases per line)");
puts("%T       reverse complement (on a single line)");
puts("%G       reverse complement (60 bases per line)");
puts("%R       sequence rank (in the universe)");
puts("%L       sequence length");
puts("%A       annotation (if any)");
puts("%C       word count (useful in single-word mode)");
puts("%n       newline");
puts("%t       tab");
puts("%%       percent sign");
puts("Anything else is copied verbatim");
puts("Example: the format strings \">%I%n%S%n\" and \">%I%n%F%n\" encode FASTA format");
return SYL_DONE;
      endarg()

      uniarg("--amoixa")
puts("--fgets              use fgets rather than fread for file read (not interesting)");
      endarg()
      uniarg("--fgets")          unset(MODE_FREAD);            endarg()

      uniargs("-h", "--apropos")
puts("--version            version information, acknowledgements");
puts("-fasta <fname>       fasta file (- for stdin)");
puts("-o <fname>           result file (- for stdout (default))");
puts("-k <int>             word length K (default 6)");
puts("-m <int>             conditionalize on word length <int>");
puts("-tt <int>            output <int> most over/under-represented words");
puts("");

puts("-subset <fname>      file containing subset IDs. Occurrence = rank if applicable");
puts("-universe <fname>    same as above, disregard anything else");
puts("-grow <int>          analyze nested leading bins with size increment <int>");
puts("");

puts("-words <fname>       read words to search from file");
puts("-w w1[,w2[,w3[..]]]  comma-separated list of words to include");
puts("-xw w1[,w2[,w3[..]]] comma-separated list of words to exclude");
puts("-entropy-gq <num>    filter out words below entroy");
puts("-entropy <word>      print entropy of word and exit");
puts("--bag                add bags for all selected words");
puts("");

puts("-step <int>          analyze consecutive bins of same size <int>");
puts("-step-times <int>    analyze <int> consecutive bins of same size");
puts("-grow-times <int>    analyze <int> nested bins of growing size");
puts("-shift <int>         ignore the first and last <int> sequences (from query)");
puts("-lshift <int>        ignore the first <int> sequences (from query)");
puts("-rshift <int>        ignore the last <int> sequences (from query)");
puts("-do <int>            stop after computing bin <int>");
puts("--reverse            reverse the input order (query file or universe)");
puts("--sort-by-length     sort sequences by length (cf. --reverse)");
puts("--bin-by-seqs        use sequence as binning unit rather than nucleotides");
puts("");

puts("-write-missing <fname> write missing IDs to file");
puts("-write-found <fname> write found IDs to file");
puts("-write-query <fname> write deduplicated IDs to file");
puts("");

puts("--over               only do overrepresentation (right/upper tail)");
puts("--under              only do underrepresentation (left/lower tail)");
puts("--logfold            compute logfold over/under-representation");
puts("--none               skip representation analysis entirely");
puts("");

puts("--print-id           append, for each bin, concatenation of all (gene) IDs");
puts("--two-strands        consider both strands (doubles universe/subset)");
puts("--no-binomial        no binomial");
puts("-trial-extremes <int> run repeated trials, dump max bin per word, <int> times");
puts("-trial-all <int>     run repeated trials, all words, all bins, <int> times");
puts("-co <num>            p-value filter");
puts("-oob <num>           cap log-transformed P-values at <num>");
puts("");

puts("-dump-expected <prefix> prefix for dumping per-bin expected frequencies");
puts("-u <num>             fake universe size when computing expected values");
puts("-u-times <num>       fake universe size as multiple of subset size");
puts("-fake-unit-size <num>fake unit size / simulate equal sequence lengths");
puts("--no-anchor          turn off expected count modulation (relevant to -m)");
puts("");

puts(" ^_^ in dyad mode use even -k to indicate summed dyad length");
puts("--dyad-repeat        require words to be a 2-repeat of the same dyad");
puts("--dyad-invert        require words to be a dyad with its reverse complement");
puts("--dyad-paired        require words to be a repeat or invert dyad as above");
puts("--monad-uniform      use uniform placement model for background correction");
puts(" ^_^ employ unrestricted dyad mode by using none of the options above");
puts("     and employing at a minimum -gap-min and -gap-max");
puts("-gap-min             minimum gap size for dyad analysis");
puts("-gap-max             maximum gap size for dyad analysis");
puts("-gap-gap             minimum allowed distance between dyad hits");
puts("");

puts("-v <int>             set verbosity level");
puts("-d <int>             set debug level");
puts("-log <fname>         log file");
puts("");

puts("--ck                 count sylamer space in units of K");
puts("-length-co <int>     cap sequences at this size; use negative to cap from end");
puts("--utr-attributes     include UTR attributes {utrlen, nmasked, gccontent} in output");
puts("");

puts("-h                   this help");
puts("--apropos            this help");
puts("--stop-bg            stop after computing background");
puts("--stop-read          stop after reading fasta file");
puts("--stop-words         stop after constructing word list (and dumping it)");
puts("--long-listing       no table output, long listing of counts/scores instead");
puts("--alien-ok           do not warn for bases not in [ACGTUNX]");
puts("-read-expected <fname>  read in expected frequencies of occurrence");
puts("-threshold <num>     count per-sequence occurrences exceeding <num> (research)");
puts("--mww                do Mann-Whitney-Wilcoxon test (research)");
puts("--shuffle            shuffle universe (trial/test mode)");
puts("--clean-up           release all memory when done (unit test -- memory)");
puts("");

puts("REPEAT MASKING OPTIONS");
puts("--r2-off             do not disregard 1-shift or 2-shift repeats");
puts("-rm <int>            reduce repeats of length up to <int> (keep ulterior ends)");
puts("-rm-char <char>      repeat mask with character <char>");
puts("--rm-maskall         mask entire repeat");
puts("export SYL_FN_REPEAT=fname && sylamer [options] # write repeats to file");
puts("");

puts("COUNT MODE options (-w/-words can be useful)");
puts("--word-count         output sequence-word counts in table format");
puts("");

puts("SEQUENCE MODE options");
puts("--format-fasta       output full fasta sequence format");
puts("--format-id          output IDs only");
puts("-format <FORMAT>     output sequences as in <FORMAT> (see --format-syntax)");
puts("--format-syntax      describe -format syntax");
puts("");

puts("FILTER MODE options (-w/-words can be useful)");
puts("-filter-lq <num>     output IDs with at most <num> counts (AND over words)");
puts("-filter-gq <num>     output IDs with at least <num> counts (OR over words)");
puts("--filter-and         require all conditions be met rather than at least one");

puts("");
puts("MISC (-w/-words can be useful)");
puts("--stat               output base count and mask count");
puts("--word-stat          output expected and observed frequencies");
puts("--word-utr-stat      probably beautiful");

puts("");
puts("For any file name, use - to get stdin or stdout");
#if HAVE_ZLIB
puts("This sylamer can read gzip-compressed files");
#endif
         return SYL_DONE;
      endarg()

      optarg("-test")            par->test = thearg();                  endarg()
      optarg("-d")               debug_G = atoi(thearg());              endarg()
      optarg("-v")               verbose_G = atoi(thearg());            endarg()
      optarg("-nalloc")          n_alloc_exit_G = atoi(thearg());       endarg()
      optarg("-write-missing")   par->fnmissing = thearg();             endarg()
      optarg("-write-found")     par->fnfound = thearg();               endarg()
      optarg("-write-query")     par->fnwquery = thearg();              endarg()
      optarg("-oob")             out_of_bounds_g = atoi(thearg());      endarg()
      optarg("-tt")              par->toptable = atoi(thearg());        endarg()
#if 0
      optarg("-aa")              par->accelerate_offset = atoi(thearg()); endarg()
      optarg("-a")               par->accelerate_num = atoi(thearg());  endarg()
#endif
      optarg("-u")               par->n_universe_fake = atoi(thearg()); endarg()
      optarg("-fake-unit-size")  par->unit_size_fake = atoi(thearg());  endarg()
      optarg("-u-times")         par->utimes_fake = atoi(thearg());     endarg()
      optarg("-cap")             par->cap = atoi(thearg());             endarg()
      optarg("-length-co")       par->length_co = atoi(thearg());       endarg()
      optarg("-rm")              par->repeat_mask = atoi(thearg());     endarg()
      optarg("-threshold")       par->per_sequence = atoi(thearg());    endarg()
      optarg("-rshift")          par->shift_right = atoi(thearg());     endarg()
      optarg("-lshift")          par->shift_left = atoi(thearg());      endarg()
      optarg("-gap-min")         par->gap_min = atoi(thearg());         endarg()
      optarg("-gap-max")         par->gap_max = atoi(thearg());         endarg()
      optarg("-gap-gap")         par->gap_gap = atoi(thearg());         endarg()
      optarg("-rm-char")         rmchar_G = thearg()[0];                endarg()
      optarg("-shift")           par->shift_left = atoi(thearg()); par->shift_right = par->shift_left; endarg()
      optarg("-fasta")           par->fpfasta = myfopen(thearg(), "r"); endarg()      /* checked later */
      optarg("-log")             par->iffplog = myfopen(thearg(), "w"); endarg()      /* subsequent use is conditional */
      optarg("-dump-expected")   par->fnexp = thearg();                 endarg()
      optarg("-co")              par->cutoff = atof(thearg());          endarg()
      optarg("-w")               par->patlist = thearg();               endarg()
      optarg("-xw")              par->xpatlist = thearg();              endarg()
      optarg("-do")              par->lastbin = atoi(thearg());         endarg()
      optarg("-step-times")      par->splitsize = atoi(thearg()); set2(Mode_STEPPING); endarg()
      optarg("-step")            par->binsize = atoi(thearg()); set2(Mode_STEPPING); endarg()
      optarg("-grow")            par->binsize = atoi(thearg());         endarg()
      optarg("-grow-times")      par->splitsize = atoi(thearg());       endarg()
      optarg("-k")               par->K_g = atoi(thearg());             endarg()
      optarg("-m")               par->M_g = atoi(thearg());             endarg()
      optarg("-o")               par->fnresult = thearg();              endarg()
      optarg("-trial-all")       par->n_trial = 1 + atoi(thearg()); set(MODE_SHUFFLE); endarg()
      optarg("-trial-extremes")  par->n_trial = 1 + atoi(thearg()); main_set(MAIN_EXTREME); if (par->n_trial > 1) set(MODE_SHUFFLE); endarg()
      optarg("-filter-lq")       par->unit_at_most = atoi(thearg());  main_set(MAIN_FILTER); endarg()
      optarg("-filter-gq")       par->unit_at_least = atoi(thearg()); main_set(MAIN_FILTER); endarg()
      optarg("-entropy-gq")      par->entropy = atof(thearg());         endarg()

         optarg("-entropy")
            fprintf(stdout, "%.6f\n", get_string_entropy(thearg()));
            return SYL_DONE;
         endarg()

         optarg("-hm")  /* undocumented */
            par->hash_commit_min_g = atoi(thearg());
            if (par->hash_commit_min_g > 1 << 14) par->hash_commit_min_g =  1 << 14;
         endarg()

         optarg("-hc")  /* undocumented */
            par->hash_commit_fac_g = atoi(thearg());
            if (par->hash_commit_fac_g < 2) par->hash_commit_fac_g = 2;
         endarg()

      optarg("-read-expected")
         if (par->fnwords) {
            return syldone(SYL_FAIL, "use only one of -words/-read-expected");
         }
         par->fnwords = thearg();
         set(MODE_FREQUENCIES);
      endarg()

      optarg("-words")
         if (par->fnwords) {
            return syldone(SYL_FAIL, "use only one of -words/-read-expected");
         }
         par->fnwords = thearg();
      endarg()

      optarg("-subset")
         if (par->fpquery) {
            return syldone(SYL_FAIL, "use only one of -subset/-universe");
         }
         if (!(par->fpquery = myfopen(thearg(), "r"))) {
            return syldone(SYL_FAIL, NULL);
         }
      endarg()

      optarg("-universe")
         if (par->fpquery) {
            return syldone(SYL_FAIL, "use only one of -subset/-universe");
         }
         if (!(par->fpquery = myfopen(thearg(), "r"))) {
            return syldone(SYL_FAIL, NULL);
         }
         set(MODE_UNIVERSE);
      endarg()

      failarg()
      arg_done()  

      return SYL_SUCCESS;
   }


