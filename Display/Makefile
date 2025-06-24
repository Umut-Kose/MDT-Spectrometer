# Compiler and flags
CXX         := g++
CXXFLAGS    := -std=c++17 -g -O0 -Wall -Wextra
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --glibs) -lEve -lRGL -lGeom -lGed -lGLEW

# Targets
TARGET      := eve_display
DICT        := MyDisplayDict.cxx
HEADERS     := MyDisplay.h
SOURCES     := MyDisplay.cpp main.cpp
OBJECTS     := $(SOURCES:.cpp=.o) $(DICT:.cxx=.o)

# Rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

MyDisplayDict.cxx: $(HEADERS) LinkDef.h
	rootcling -f $@ -c $(HEADERS) LinkDef.h

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

MyDisplayDict.o: MyDisplayDict.cxx
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS) $(DICT) *.pcm *.d *~ core

distclean: clean
	rm -rf eve_display.dSYM

.PHONY: all clean distclean
