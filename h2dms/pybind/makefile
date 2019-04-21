NAME=h2hp
CXX=g++
OFLAGS=-O3
CFLAGS=$(OFLAGS) -shared -fPIC
OBJS=h2hp.o h2ham.o j2pade.o prodwf.o pybind.o
LIB=$(NAME).so

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $<

all: $(OBJS)
	$(CXX) $(CFLAGS) -o $(LIB) $(OBJS)

clean:
	rm $(OBJS) $(LIB); find -type l -delete
