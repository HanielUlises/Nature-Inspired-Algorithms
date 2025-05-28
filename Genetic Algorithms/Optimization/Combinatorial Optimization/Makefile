CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

OBJS = tsp_solver.o main.o

all: tsp_solver

tsp_solver: $(OBJS)
	$(CXX) $(CXXFLAGS) -o tsp_solver $(OBJS)

tsp_solver.o: tsp_solver.cpp tsp_solver.h
	$(CXX) $(CXXFLAGS) -c tsp_solver.cpp

main.o: main.cpp tsp_solver.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f *.o tsp_solver
