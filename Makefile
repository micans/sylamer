

GSLPREFIX = $(HOME)/local

D = -DHAVE_ZLIB=1
L = -lgsl -lgslcblas -lm -lz

##  if zlib is not on your system, outcomment the lines below
# D = 
# L = -lgsl -lgslcblas -lm

CC = gcc
C = -Wall -pedantic -g -O2



.PHONY: clean wc
   
modules_o = table.o util.o sylio.o query.o words.o unit.o array.o acc.o acciter.o interval.o interface.o
modules_c = $(modules_o:.o=.c) 
modules_h = $(modules_o:.o=.h) 

sylamer: sylamer.o $(modules_o) version.h
	$(CC) -L$(GSLPREFIX)/lib $(C) $(D) -o sylamer sylamer.o $(modules_o) $L

install:
	cp sylamer $(HOME)/local/bin

hyperdrive: hyperdrive.o
	$(CC) -L$(GSLPREFIX)/lib $(C) -o hyperdrive hyperdrive.o $L

testset: testset.o
	$(CC) -L$(GSLPREFIX)/lib $(C) -o testset testset.o $L

static: sylamer.o $(modules_o)
	$(CC) $(C) $(D) -o sylamer-static sylamer.o $(modules_o) ${GSLPREFIX}/lib/libgsl.a ${GSLPREFIX}/lib/libgslcblas.a /usr/lib64/libz.a -lm

hyperdrive.o: hyperdrive.c
	$(CC) -c hyperdrive.c -I$(GSLPREFIX)/include $(D) $(C)

testset.o: testset.c
	$(CC) -c testset.c -I$(GSLPREFIX)/include $(D) $(C)

sylamer.o: sylamer.c
	$(CC) -c sylamer.c -I$(GSLPREFIX)/include $(D) $(C)

table.o: table.c util.h
	$(CC) -c -o table.o table.c $(C)

util.o: util.c
	$(CC) -c -o util.o util.c $(C)

sylio.o: sylio.c util.h
	$(CC) -c -o sylio.o sylio.c $(D) $(C)

query.o: query.c util.h sylio.h
	$(CC) -c -o query.o query.c $(C)

words.o: words.c util.h sylio.h
	$(CC) -c -o words.o words.c $(C)

unit.o: unit.c query.h util.h
	$(CC) -c -o unit.o unit.c $(C)

array.o: array.c util.h
	$(CC) -c -o array.o array.c $(C)

acc.o: acc.c util.h array.h unit.h
	$(CC) -c -o acc.o acc.c $(C)

acciter.o: acciter.c util.h acc.h array.h
	$(CC) -c -o acciter.o acciter.c $(C)

interval.o: interval.c util.h
	$(CC) -c -o interval.o interval.c $(C)

interface.o: interface.c util.h unit.h array.h acc.h version.h
	$(CC) -c -o interface.o interface.c $(C) $(D)

clean:
	rm -f sylamer sylamer.o $(modules_o)

