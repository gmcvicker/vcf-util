# Edit the following variables as needed
CC=gcc
LIB=-lz

INCLUDE=
CFLAGS=-g $(INCLUDE)

objects=vcf.o util.o memutil.o err.o chrom.o

default: all

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDE)

vcfmerge: $(objects) vcfmerge.c
	$(CC) $(CFLAGS) -o $@ $(objects) vcfmerge.c $(LIBSHDF) $(LIB)

all:  $(objects) vcfmerge

clean:
	rm -f $(objects) vcfmerge
