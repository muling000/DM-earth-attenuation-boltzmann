CXX = g++
CXXFLAGS = -I../.. $(shell pkg-config --cflags spdlog) -Wall -O2 -std=c++17
LDFLAGS = $(shell pkg-config --libs spdlog)

run.exe : trajectory.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -rf *.exe

.PHONY:
	clean
