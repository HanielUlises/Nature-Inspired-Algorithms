TARGET = genetic_algorithm

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

SRCS = genetic_algo.cpp PM/polynomial_mut.cpp SBX/sbx.cpp ../../Utils/ObjectiveFunctions.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)

run: all
	./$(TARGET)