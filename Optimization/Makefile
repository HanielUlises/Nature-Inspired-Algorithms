# Compiler
CXX = clang++

# Compiler flags
CXXFLAGS = -std=c++20 -Wall

# Target executable name
TARGET = program.exe

# Source files
SOURCES = main.cpp

# Object files (automatically deduced from source files)
OBJECTS = $(SOURCES:.cpp=.o)

# Header files
HEADERS = function_test_1.h

# Default rule to build the executable
all: $(TARGET)

# Rule to link the program
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

# Rule to compile the source files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(TARGET) $(OBJECTS)

# Prevent make from doing something with a file named all or clean
.PHONY: all clean
