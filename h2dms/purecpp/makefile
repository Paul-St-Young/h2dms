NAME=mmh
CC=g++
OBJS=main.o h2hp.o
#OFLAGS=-O0 -g
OFLAGS=-O3
CFLAGS=-c $(OFLAGS) -I. -std=c++11

%.o : %.cpp
	$(CC) $(CFLAGS) $<

all: $(OBJS)
	$(CC) -o exec $(OFLAGS) $(OBJS)
clean:
	rm *.o exec
