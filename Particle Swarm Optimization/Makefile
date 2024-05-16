CXX = clang++
CXXFLAGS = -std=c++20 -Wall -I"../Utils"
LDFLAGS =

TARGET = PSO.exe

SRC_DIR = .
OBJ_DIR = obj
UTILS_DIR = ../Utils

SRCS = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(UTILS_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cpp)) \
       $(patsubst $(UTILS_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(wildcard $(UTILS_DIR)/*.cpp))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $<"
	@if not exist "$(OBJ_DIR)" mkdir "$(OBJ_DIR)"
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(UTILS_DIR)/%.cpp
	@echo "Compiling $<"
	@if not exist "$(OBJ_DIR)" mkdir "$(OBJ_DIR)"
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	del /Q /F $(OBJ_DIR)\\* $(TARGET)

.PHONY: all clean
