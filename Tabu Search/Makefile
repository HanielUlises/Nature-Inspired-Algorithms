CXX = c++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -g

TARGET = main.exe
SOURCES = main.cpp tabu_search.cpp

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

clean:
	rm -f $(TARGET)
