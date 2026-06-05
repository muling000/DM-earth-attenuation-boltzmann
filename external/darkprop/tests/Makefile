CXX ?= g++

GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LDFLAGS := $(shell gsl-config --libs)

CXXFLAGS = -Wall -std=c++17 -g $(GSL_CFLAGS) $(shell pkg-config --cflags spdlog)

LDFLAGS = $(GSL_LDFLAGS) $(shell pkg-config --libs spdlog)

HEADERS = $(wildcard *.hpp ../darkprop/*.hpp)
SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%.o,${SRCS})


test : $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
	./test

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f test *.o

.PHONY:
	clean
