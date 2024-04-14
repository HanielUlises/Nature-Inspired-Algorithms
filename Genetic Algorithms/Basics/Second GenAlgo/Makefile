CXX = g++
CXXFLAGS = -std=c++20 -I/usr/include/python3.10 -I/home/haniel/Libraries -I/home/haniel/.local/lib/python3.10/site-packages/numpy/core/include
LDFLAGS = -lpython3.10

TARGET = main

SRC = main.cpp EvolutionaryAlgorithm.cpp ploteado.cpp

OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
