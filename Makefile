CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
LDFLAGS = -lm

# Source files
SOURCES = Nanopyramid.cpp metropolis.cpp coordinates.cpp geometry.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = Nanopyramid

# Default target
all: $(EXECUTABLE)

# Link the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files
%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Special case for main file if it doesn't have a header
Nanopyramid.o: Nanopyramid.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

# Rebuild everything
rebuild: clean all

.PHONY: all clean rebuild
