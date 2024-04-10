all: magic_square

magic_square: main.o EvolutionaryAlgorithm.o
	g++ -o magic_square main.o EvolutionaryAlgorithm.o

main.o: main.cpp
	g++ -c main.cpp

EvolutionaryAlgorithm.o: EvolutionaryAlgorithm.cpp
	g++ -c EvolutionaryAlgorithm.cpp

clean:
	rm *.o magic_square
