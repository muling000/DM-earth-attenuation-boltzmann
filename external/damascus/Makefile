
#Compiler and compiler flags
CXX := mpic++
CXXFLAGS := -Wall -std=c++14
EIGEN_CFLAGS ?= $(shell pkg-config --cflags eigen3 2>/dev/null)
ifeq ($(strip $(EIGEN_CFLAGS)),)
EIGEN_CFLAGS := -I/usr/local/include/eigen3
endif
LIBCONFIG_LIBS ?= $(shell pkg-config --libs libconfig++ 2>/dev/null)
ifeq ($(strip $(LIBCONFIG_LIBS)),)
LIBCONFIG_LIBS := -lconfig++
endif
LIBCONFIG_CFLAGS ?= $(shell pkg-config --cflags libconfig++ 2>/dev/null)
LIB := $(LIBCONFIG_LIBS)
INC := -I include $(EIGEN_CFLAGS) $(LIBCONFIG_CFLAGS)

#Directories
SRCDIR := src
BUILDDIR := build
TARGETDIR := bin

#Targets
TARGETS = $(TARGETDIR)/DaMaSCUS-Simulator $(TARGETDIR)/DaMaSCUS-Analyzer

#Source files
SRCEXT := cpp
COMMONSRC := $(shell find $(SRCDIR) -maxdepth 1 -type f -name '*.$(SRCEXT)')
SIMSRC :=$(COMMONSRC) $(shell find $(SRCDIR)/simulation -type f -name *.$(SRCEXT))
ANASRC :=$(COMMONSRC) $(shell find $(SRCDIR)/analysis -type f -name *.$(SRCEXT))

#Object files
SIMOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SIMSRC:.$(SRCEXT)=.o))
ANAOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(ANASRC:.$(SRCEXT)=.o))

.PHONY: all simulator analyzer clean codecov


all: CXXFLAGS += -O2 
all: $(TARGETS)

test: LIB+= --coverage
test: $(TARGETS)

simulator: CXXFLAGS += -O2
simulator: $(TARGETDIR)/DaMaSCUS-Simulator

analyzer: CXXFLAGS += -O2
analyzer: $(TARGETDIR)/DaMaSCUS-Analyzer

$(TARGETDIR)/DaMaSCUS-Simulator: $(SIMOBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)

$(TARGETDIR)/DaMaSCUS-Analyzer: $(ANAOBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)/simulation/
	@mkdir -p $(BUILDDIR)/analysis/
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

codecov:
	

clean:
	find . -type f -name '*.o' -delete
	find . -type f -name '*.gcno' -delete
	find . -type f -name '*.gcda' -delete
	find . -type f -name '*.gcov' -delete
	rm -f $(TARGETS)
