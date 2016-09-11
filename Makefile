# Makefile for SWEMBL
CC = gcc
CFLAGS = -g -Wall -O3
LFLAGS = -lz -lm
OBJS = main.o IO.o calc.o stack.o summit.o refcalc.o wiggle.o overlap.o
EFILE = SWEMBL

default: ${EFILE}

%.o: %.c SWEMBL.h
	$(CC) $(CFLAGS) -c $<

$(EFILE): $(OBJS)
	$(CC) $(CFLAGS) -o $(EFILE) $(OBJS) $(LFLAGS) 

clean:
	rm *.o
