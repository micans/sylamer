

#include "sylio.h"
#include "util.h"
#include "inttypes.h"
#include "acc.h"

#include <string.h>
#include <ctype.h>

#ifndef DEBUG
#define DEBUG 0
#endif

   /*
   */
const char* mygets(INPUT fp, dim* length)
#define BUFFER_MYGETS_SIZE 4097
   {  static char mygets_buf[BUFFER_MYGETS_SIZE]
   ;  static dim  mygets_ofs = 0
   ;  static dim  mygets_n_read = 0

   ;  dim ofs_ret
   ;  size_t n_spn = 0

   ;  if (mygets_ofs >= mygets_n_read)
      {
#if HAVE_ZLIB
         {  ofs n_read = gzread(fp, mygets_buf, BUFFER_MYGETS_SIZE-1);
         ;  if (n_read <=0)
            return NULL
         ;  mygets_n_read = n_read
      ;  }
#else
         if
         (! (  mygets_n_read
            =  fread(mygets_buf, 1, BUFFER_MYGETS_SIZE-1, fp)
         )  )
         return NULL
#endif
      ;  mygets_ofs = 0
      ;  mygets_buf[mygets_n_read] = '\0'   /* use this as sentinel below */
   ;  }

      {  char* o = mygets_buf + mygets_ofs
      ;  while (o[0] && !(8 & basemap[(unsigned char) o[0]]))
         o++
                        /* v caller is responsible for removing trailing whitespace */
      ;  if ((unsigned char) o[0] == '\r' && (unsigned char) o[1] == '\n')
         o++
      ;  n_spn = (o - mygets_buf) - mygets_ofs
   ;  }

            /* If \r or \n was found convert to \n always.  We may splice a
             * DOS file as [...\r][\n....].  If the condition below is not met
             * it means none of the retrieved characters was a newline.
            */
   ;  if (mygets_n_read && n_spn + mygets_ofs < mygets_n_read)
      {  mygets_buf[mygets_ofs + n_spn] = '\n'
      ;  *length = n_spn + 1
   ;  }
      else
      *length = n_spn

   ;  ofs_ret = mygets_ofs
   ;  mygets_ofs += length[0]

   ;  return mygets_buf+ofs_ret
;  }


/* fixme - EOF case */


   /* This buffers read with fread. It is quite fast, but the gain over fgets
    * is not big compared with other costs we have.
    * Returned line is stripped of newline.
    *
    * If instructed to flush,
    *    return NULL and set status to SYL_IGNORE when sequence is found.
    *
    * Otherwise, when it finds an identifier, 
    *    it will return line and set SYL_SUCCESS.
    *
    * At EOF it sets SYL_DONE.
   */
char* fread_line_dos_mac_unx
(  INPUT fp
,  int*  have_id
,  int   flush_sequence
,  dim*  ll
,  sylstatus* status
)
   {  char dummy_line[] = "thedishrawnawaywiththespoon"
   ;  char* line = NULL, *line_realloc = dummy_line   /* bit hackish, but needed */
   ;  dim   linelen = 0, n_chunk = 0
   ;  int   delta = 0, eol = 0
   ;  *ll = 0

   ;  *status = SYL_SUCCESS

   ;  if (have_id)
      *have_id = 0

   ;  while (!eol)
      {  const char* coffset = 0
      ;  dim length = 0, addlength = 0

      ;  if (!(coffset = mygets(fp, &length)))
         {        /* nothing was read .... */
#if HAVE_ZLIB
            if (!gzeof(fp))
#else
            if (!feof(fp))
#endif
            status[0] = arrr("line input error")
         ;  else
            status[0] = SYL_DONE
         ;  break
      ;  }

                  /* We avoid copying and mallocing when instructed to flush.
                   * However, flushing is aborted when a new identifier is foud.
                   * So we turn off flushing if we have '>'.
                  */
         if (!n_chunk && have_id && (unsigned char) coffset[0] == '>')
         {  delta = 1
         ;  *have_id = 1
         ;  flush_sequence = 0
      ;  }
         else delta = 0

      ;  n_chunk++
      ;  if ((unsigned char) coffset[length-1] == '\n')
         {  eol = 1
                  /* Below mainly is to get rid of DOS \r in \r\n endings,
                   * but will also remove trailing spaces.
                   * This branch is not triggered in the EOF case,
                   * but that is where we draw the line (noteme/fixme/fixmenot).
                  */
         ;  while(length > 1 && isspace((unsigned char) coffset[length-2]))
            length--
      ;  }

                  /* When flushing, all reallocing and length increases
                   * are avoided. We will return NULL and SYL_IGNORE.
                  */
         if (flush_sequence)
         continue

      ;  addlength = length - delta - eol       /* this excludes space for '\0' */
      ;  if (!(line_realloc = myrealloc(line, linelen + addlength + 1)))
         {  status[0] = SYL_FAIL
         ;  break
      ;  }
         line = line_realloc
      ;  memcpy(line+linelen, coffset+delta, addlength)
      ;  line[linelen+addlength] = '\0'
      ;  linelen += addlength
   ;  }

      if (status[0] == SYL_FAIL)
      {  if (line)
         free(line)
      ;  line = NULL
      ;  linelen = 0
   ;  }

      if (!status[0] && flush_sequence)
      *status = SYL_IGNORE

   ;  *ll = linelen
   ;  return line
;  }




   /* This is a wrapper. As we cater for all sorts of line endings,
    * we might end up with empty lines that are not really empty lines.  We
    * just accept any empty line and throw it away.  It is too much headache to
    * sort out 'real' empty lines from \r\n\r\r\n, (<--\r][\n-->) splice events
    * and related artifacts. Some of this headache is due to how I've setup the
    * (buffered) reading routines.  Refer also to implementation notes.
    * The flushing modality makes the code even less legible, unfortunately. fixme/cleanup.
   */
char* fread_line
(  INPUT fp
,  int* have_id
,  int flush_sequence
,  dim* linelen
,  dim* linecount
,  sylstatus * status
)
   {  char* line
                                 /* hierverder: check below gives headaches */
   ;  do
      {  line = fread_line_dos_mac_unx(fp, have_id, flush_sequence, linelen, status)
      ;  linecount[0]++
                                 /* if have_id && have_id[0] return empty identifier */
      ;  if (line && !strlen(line) && (!have_id || !*have_id))
         free(line)
      ;  else
         break
   ;  }
      while (1)
   ;  return line
;  }


char* fgets_line(INPUT fp, int* have_id, dim* ll, sylstatus* status)
#define BUFFER_FGETS_SIZE 4096
   {  static char rline_buf[BUFFER_FGETS_SIZE]
   ;  char dummy[] = "dummy"
   ;  char* line = NULL, *line_realloc = dummy
   ;  dim linelen = 0, buflen = 0, n_chunk = 0, addlen = 0
   ;  int delta = 0, eol = 0

   ;  if (have_id)
      *have_id = 0

   ;  *ll = 0

   ;  while (!eol)
      {
#if HAVE_ZLIB
         if (!gzgets(fp, rline_buf, BUFFER_FGETS_SIZE))
         {  if (!gzeof(fp))
            *status = arrr("gzgets error")
         ;  break
      ;  }
#else         
         if (!fgets(rline_buf, BUFFER_FGETS_SIZE, fp))
         {  if (!feof(fp))
            *status = arrr("fgets error")
         ;  break
      ;  }
#endif

         if (!n_chunk && have_id && (unsigned char) rline_buf[0] == '>')
         {  delta = 1
         ;  *have_id = 1
      ;  }
         else
         delta = 0

      ;  buflen = strlen(rline_buf)
      ;  if ((unsigned char) rline_buf[buflen-1] == '\n')
         eol = 1

      ;  addlen = buflen-delta-eol
      ;  line_realloc = myrealloc(line, linelen + addlen + 1)
      ;  if (!line_realloc)
         break
      ;  line = line_realloc

      ;  memcpy(line + linelen, rline_buf+delta, addlen)
      ;  line[linelen + addlen] = '\0'
      ;  linelen += addlen
      ;  n_chunk++
   ;  }
#if DEBUG & DEBUG_READ
miaow("read line [%s]", line);
#endif

      if (*status || !line_realloc)
      {  if (line)
         free(line)
      ;  *status = SYL_FAIL
      ;  line = NULL
      ;  linelen = 0
   ;  }
      *ll = linelen
   ;  return line
;  }


   /* This returns void* because it may, for input streams either return a
    * FILE* or a gzFile type, depending on HAVE_ZLIB.  This type is defined as
    * INPUT*.  However, we cannot return INPUT* as for output streams myfopen
    * always returns a FILE*.  HAVE_ZLIB makes things slightly more difficult
    * to read.
   */
void* myfopen (const char* fname, const char* mode)
   {  FILE* fp = NULL
   ;  if (!strcmp(fname, "-"))
      {  if (!strcmp(mode, "w"))
         fp = stdout
      ;  else if (!strcmp(mode, "r"))
         fp = stdin
      ;  else
         miaow("fire: unrecognized open mode <%s>", mode)     /* fp now NULL */
   ;  }
      else if (!strcmp(fname, "stderr") && !strcmp(mode, "w"))
      fp = stderr
   ;  else
      {  fp = fopen(fname, mode)
      ;  if (!fp)
         miaow("cannot open file <%s> in mode <%s>", fname, mode)
   ;  }

#if HAVE_ZLIB
      if (fp && !strcmp(mode, "r"))
      {  int fn = fileno(fp)
      ;  if (fn < 0)
         miaow("bad fileno (file %s)", fname)
      ;  return fn < 0 ? NULL : gzdopen(fn, mode)
   ;  }
#endif
      return fp
;  }


void myfzclose(INPUT* fpp)
   {  if(*fpp)
      {
#if HAVE_ZLIB
         gzclose(*fpp)
#else
         fclose(*fpp)
#endif
   ;  }
      *fpp = NULL
;  }

void myfclose(FILE** fpp)
   {  if(*fpp)
      fclose(*fpp)
   ;  *fpp = NULL
;  }


