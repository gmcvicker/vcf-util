# Edit the following variables as needed
CC=gcc
LIB=-lz

INCLUDE=
CFLAGS=-g $(INCLUDE)

objects=vcf.o util.o memutil.o err.o

default: all

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDE)

all:  $(objects)

clean:
	rm -f $(objects)
