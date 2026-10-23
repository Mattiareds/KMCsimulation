CXX = g++
CXXFLAGS = -std=c++17 -Wall
LDFLAGS = -lm

# make debug=1 for activate
ifdef debug
    CXXFLAGS += -g3 -O2 -fsanitize=address -fsanitize=undefined
    LDFLAGS += -fsanitize=address -fsanitize=undefined
else
    CXXFLAGS += -O2
endif

SOURCES = Nanopyramid.cpp metropolis.cpp coordinates.cpp geometry.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = Nanopyramid

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
