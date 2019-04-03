
#ifndef  include_words
#define  include_words

#include "inttypes.h"
#include "sylio.h"

typedef struct
{  unsigned*   words
;  dim         n_words
;  dim         n_words_alloc
;
}  wlist       ;

sylstatus get_wordlist
(  INPUT fpwords
,  wlist*  wl
,  float** oracle_pp
,  unsigned* K_pp          /* check for consistency if *K_pp != 0, otherwise set it */
)  ;

sylstatus get_dyadlist
(  unsigned modes
,  wlist* wl
,  unsigned Kh             /* dyad length (half of K). */
)  ;

sylstatus get_bypattern
(  const char* list
,  wlist* wl
,  unsigned* K_found
)  ;

sylstatus get_bags
(  wlist* wl
,  unsigned k
)  ;


void        wlist_free(wlist* wl);
sylstatus   wlist_full(wlist* wl, unsigned k);
void        wlist_remove(wlist* wl, wlist* wldelete);
dim         wlist_uniq(wlist* wl);
dim         wlist_entropy(wlist* wl, unsigned K, double e);

double get_string_entropy (const char* word);
double get_word_entropy(unsigned w, unsigned K);

#endif

