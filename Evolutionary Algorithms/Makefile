CXX = clang++
CXXFLAGS = -std=c++20

TARGET = main.exe

SRC = main.cpp Evolutionary_Strategy.cpp Differential_Evolution.cpp ObjectiveFunctions.cpp

OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
