#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "Genetic_Algorithm.h"
#include "Solution.h"

using namespace std;

struct City {
    double x, y;
    City(double x_, double y_) : x(x_), y(y_) {}
};

double calculateDistance(const City& c1, const City& c2) {
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return sqrt(dx * dx + dy * dy);
}

double calculateTotalDistance(const vector<int>& route, const vector<City>& cities) {
    double total = 0;
    for (size_t i = 0; i < route.size() - 1; i++) {
        total += calculateDistance(cities[route[i]], cities[route[i + 1]]);
    }
    total += calculateDistance(cities[route.back()], cities[route[0]]);
    return total;
}

int main() {
    srand(static_cast<unsigned>(time(0)));
    const int NUM_CITIES = 10;
    vector<City> cities;
    for (int i = 0; i < NUM_CITIES; i++) {
        cities.emplace_back(rand() % 100, rand() % 100);
    }
    const int POPULATION_SIZE = 100;
    const int GENERATIONS = 500;
    const double MUTATION_RATE = 0.01;
    Genetic_Algorithm ga(POPULATION_SIZE, NUM_CITIES, MUTATION_RATE);
    cout << "Starting TSP optimization..." << endl;
    for (int gen = 0; gen < GENERATIONS; gen++) {
        vector<Solution>& population = ga.getPopulation();
        for (Solution& sol : population) {
            double distance = calculateTotalDistance(sol.getRoute(), cities);
            sol.setFitness(1.0 / distance);
        }
        ga.evolve();
        if (gen % 50 == 0) {
            Solution best = ga.getBestSolution();
            double bestDistance = calculateTotalDistance(best.getRoute(), cities);
            cout << "Generation " << gen << ": Best Distance = " << bestDistance << endl;
        }
    }
    Solution bestSolution = ga.getBestSolution();
    vector<int> bestRoute = bestSolution.getRoute();
    double finalDistance = calculateTotalDistance(bestRoute, cities);
    cout << "\nFinal Results:" << endl;
    cout << "Best Route: ";
    for (int city : bestRoute) {
        cout << city << " ";
    }
    cout << endl;
    cout << "Total Distance: " << finalDistance << endl;
    cout << "\nCity Coordinates:" << endl;
    for (int i = 0; i < NUM_CITIES; i++) {
        cout << "City " << i << ": (" << cities[i].x << ", " << cities[i].y << ")" << endl;
    }
    return 0;
}