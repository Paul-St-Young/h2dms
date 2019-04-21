NAME=mmh
CC=g++
OBJS=main.o h2hp.o
HEADERS=common.h h2wf.h
#OFLAGS=-O0 -g
OFLAGS=-O3
CFLAGS=-c $(OFLAGS) -I. -std=c++11

%.o : %.cpp
	$(CC) $(CFLAGS) $<

%.o : %.cpp %.h
	$(CC) $(CFLAGS) $<

all: $(OBJS) $(HEADERS)
	$(CC) -o exec $(OFLAGS) $(OBJS)
clean:
	rm *.o exec
