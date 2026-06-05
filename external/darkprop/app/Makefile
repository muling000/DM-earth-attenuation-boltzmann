CXX = mpicxx

CXXFLAGS = -Wall -O2 -std=c++17 -I.. -I/usr/include/hdf5/openmpi \
		   -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi $(shell pkg-config --cflags spdlog)

LDFLAGS = -lgsl -lgslcblas -lhdf5 $(shell pkg-config --libs spdlog)

HEADERS = $(wildcard ../darkprop/*.hpp *.hpp *.h)

darkprop.exe : main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f darkprop.exe
.PHONY:
	clean
