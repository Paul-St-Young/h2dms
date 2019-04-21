NAME=mmh
CC=g++
OBJS=main.o h2hp.o h2ham.o j2pade.o prodwf.o
HEADERS=common.h h2wf.h h2hp.h h2ham.h j2pade.h prodwf.h
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
