CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic
LDFLAGS =
OBJ = GeneticAlgorithm.o main.o Solution.o Utils.o
DEPS = GeneticAlgorithm.h Solution.h Utils.h
TARGET = genAlgo

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f $(OBJ) $(TARGET)
