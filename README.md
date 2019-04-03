# Sylamer
  Fast assessment of microRNA binding and siRNA off-target effects from expression data

# Installation
  You'll need to have GSL (Gnu Scientific Library) installed on your system,
  version 1.11 or later.  Version 1.10 is the first GSL version that contains
  the required hypergeometric functionality. However, it has a bug that may
  cause negative probabilities to be computed (as a result of computing the
  wrong tail). This bug was fixed for 1.11, which is by now an ancient version.
  In the Makefile, the variable GSLPREFIX needs to be set to the proper value.
  It should be such that

  `$(GSLPREFIX)/lib`  contains the GSL libraries (libgsl.a libgslcblas.a and/or
     the .so variants of these).

  `$(GSLPREFIX)/include` contains the GSL header files -- these will be in the
     directory `$(GSLPREFIX)/include/gsl`

  If, for some reason, on your system the GSL library files and header files
  do not share the same prefix as described above, make the necessary
  modifications in the Makefile -- it should be straightforward.
  After you have done this, typing 'make' should create a functioning
  sylamer executable.

  By default the Makefile is set up such that the presence of zlib libraries
  on your system is assumed. Zlib enables Sylamer to read and write
  gzip-compressed files. The Makefile contains a brief description of what
  to do in case zlib is not on your system. Please note that it could
  be possible that zlib is on your system but for some reason not
  in the default include paths for your compiler. In that case, the
  solution is to provide the right include path and linker options
  to the compiler. Get in touch if you find yourself in this situation.

  The manual page is available in HTML (sylamer.html), PDF (sylamer.pdf),
  PostScript (sylamer.ps), text format (sylamer.txt) and as a UNIX manual page
  (troff format, sylamer.1).

     sylamer --help or sylamer -h gives a synopsis of available options.

  If you have questions or comments about this program please send a message
  to stijn at sanger dot ac dot uk .

