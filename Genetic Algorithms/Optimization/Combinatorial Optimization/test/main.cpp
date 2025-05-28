// main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>

using namespace std;

// Estructura que define cada ciudad con ventana de tiempo.
struct City {
    string name;
    double windowOpen;
    double windowClose;
};

// Estructura que representa un individuo (ruta) del AG.
struct Individual {
    vector<int> route;
    double fitness;    
};

// Datos globales: matriz de distancias y las ciudades con ventanas.
const vector<vector<double>> distances = {
    {0, 61.82, 18.54, 37.52, 54.08, 1.88, 59.98, 32.82, 69.42, 36.76, 60.26},
    {61.82, 0, 50.84, 33.62, 7.5, 59.88, 2.76, 28.84, 7.78, 28.14, 5.8},
    {18.54, 50.84, 0, 26.74, 43.38, 18.6, 49.28, 22.0, 58.7, 23.36, 49.3},
    {37.52, 33.62, 26.74, 0, 26.16, 35.56, 32.06, 4.8, 41.5, 3.26, 32.08},
    {54.08, 7.5, 43.38, 26.16, 0, 52.06, 7.32, 21.38, 15.34, 20.68, 5.92},
    {1.88, 59.88, 18.6, 35.56, 52.06, 0, 57.96, 30.86, 67.38, 34.8, 58.3},
    {59.98, 2.76, 49.28, 32.06, 7.32, 57.96, 0, 27.28, 10.62, 26.58, 6.76},
    {32.82, 28.84, 22.0, 4.8, 21.38, 30.86, 27.28, 0, 36.72, 4.02, 27.3},
    {69.42, 7.78, 58.7, 41.5, 15.34, 67.38, 10.62, 36.72, 0, 36.02, 12.14},
    {36.76, 28.14, 23.36, 3.26, 20.68, 34.8, 26.58, 4.02, 36.02, 0, 26.6},
    {60.26, 5.8, 49.3, 32.08, 5.92, 58.3, 6.76, 27.3, 12.14, 26.6, 0}
};

const vector<City> cities = {
    {"New York",      0,   numeric_limits<double>::infinity()},  
    {"Los Angeles",  50,  90},
    {"Chicago",      15,  25},
    {"Houston",      30,  55},
    {"Phoenix",      15,  75},
    {"Philadelphia",  5,  35},
    {"San Diego",   150, 200},
    {"Dallas",       25,  50},
    {"San Francisco",65, 100},
    {"Austin",      120, 150},
    {"Las Vegas",    30,  85}
};

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// Función que evalúa la aptitud (fitness) de una ruta.
// Se considera el tiempo acumulado (tiempo de viaje + esperas) y se penaliza si se llega tarde.
// Se asume que la ruta es una secuencia de índices de ciudades y se recorre en orden.
double evaluateFitness(const vector<int>& route) {
    double time = 0.0;
    double penaltyFactor = 10.0; 
    double totalCost = 0.0;
    
    // Suponemos que la ruta se inicia en la primera ciudad de la permutación.
    // Se evalúa cada llegada.
    int currentCity = route[0];
    // Si se llega antes de la ventana, se espera.
    if(time < cities[currentCity].windowOpen)
        time = cities[currentCity].windowOpen;
    // Si se llega después de la ventana, se penaliza.
    if(time > cities[currentCity].windowClose)
        totalCost += penaltyFactor * (time - cities[currentCity].windowClose);
    
    // Recorrer el resto de las ciudades.
    for (size_t i = 1; i < route.size(); i++) {
        int nextCity = route[i];
        // Sumar tiempo de viaje
        time += distances[currentCity][nextCity];
        
        // Ajuste por ventana de tiempo:
        if(time < cities[nextCity].windowOpen) {
            time = cities[nextCity].windowOpen; // esperar hasta la apertura
        }
        if(time > cities[nextCity].windowClose) {
            totalCost += penaltyFactor * (time - cities[nextCity].windowClose);
        }
        
        currentCity = nextCity;
    }
    // La función objetivo combina el tiempo total (incluyendo esperas) y las penalizaciones.
    totalCost += time;
    return totalCost;
}

// Heurística de “Remoción de Abruptos” a través de una mejora 2-opt (básica).
// Se recorre la ruta y se intercambian segmentos si se mejora la aptitud.
void removeAbruptChanges(vector<int>& route) {
    bool improved = true;
    while (improved) {
        improved = false;
        for (size_t i = 1; i < route.size() - 1; i++) {
            for (size_t j = i + 1; j < route.size(); j++) {
                vector<int> newRoute = route;
                // Realizar intercambio 2-opt: invertir el segmento entre i y j.
                reverse(newRoute.begin() + i, newRoute.begin() + j + 1);
                if (evaluateFitness(newRoute) < evaluateFitness(route)) {
                    route = newRoute;
                    improved = true;
                }
            }
        }
    }
}

// Operador de cruce por "Cycle Crossover (CX)" para permutaciones.
// Se basa en copiar en el hijo la posición de uno de los padres para las posiciones
// que pertenezcan a un ciclo determinado y, en las demás, se toma del otro padre.
vector<int> cycleCrossover(const vector<int>& parent1, const vector<int>& parent2) {
    size_t n = parent1.size();
    vector<int> child(n, -1);
    vector<bool> taken(n, false);
    
    // Inicia el ciclo desde la posición 0.
    int index = 0;
    do {
        child[index] = parent1[index];
        taken[index] = true;
        // Buscar la posición en parent1 donde se encuentra el elemento de parent2 en la posición actual.
        int element = parent2[index];
        index = find(parent1.begin(), parent1.end(), element) - parent1.begin();
    } while(!taken[index]);
    
    // Para los elementos no copiados (fuera del ciclo) se toman de parent2.
    for (size_t i = 0; i < n; i++) {
        if(child[i] == -1)
            child[i] = parent2[i];
    }
    return child;
}

// Función auxiliar para generar una permutación aleatoria (ruta) de los índices de ciudades.
vector<int> generateRandomRoute(int n) {
    vector<int> route(n);
    for (int i = 0; i < n; i++)
        route[i] = i;
    shuffle(route.begin(), route.end(), rng);
    return route;
}

int main() {
    // Parámetros del algoritmo
    int populationSize = 10;
    int generations = 100;           // Número de generaciones (iteraciones)
    double pm = 0.1;                 // Probabilidad de inyección de una ruta aleatoria (mutación global)
    int numCities = cities.size();
    
    // Generación aleatoria de la población inicial.
    vector<Individual> population;
    for (int i = 0; i < populationSize; i++) {
        Individual indiv;
        indiv.route = generateRandomRoute(numCities);
        removeAbruptChanges(indiv.route); // aplica heurística a cada ruta
        indiv.fitness = evaluateFitness(indiv.route);
        population.push_back(indiv);
    }
    
    // Ciclo principal del AG.
    for (int gen = 0; gen < generations; gen++) {
        // Para almacenar la nueva población.
        vector<Individual> newPopulation;
        
        // Mientras la nueva población no alcance el tamaño deseado
        while(newPopulation.size() < populationSize) {
            // Seleccionar dos padres de forma aleatoria.
            uniform_int_distribution<int> distIndex(0, population.size() - 1);
            const vector<int>& parentRoute1 = population[distIndex(rng)].route;
            const vector<int>& parentRoute2 = population[distIndex(rng)].route;
            
            // Aplicar ciclo de cruce (CX) para generar un descendiente.
            vector<int> childRoute = cycleCrossover(parentRoute1, parentRoute2);
            removeAbruptChanges(childRoute); // aplicar heurística al descendiente
            double childFitness = evaluateFitness(childRoute);
            
            // Crear “familia” (dos padres y el descendiente)
            vector<Individual> family;
            family.push_back({parentRoute1, evaluateFitness(parentRoute1)});
            family.push_back({parentRoute2, evaluateFitness(parentRoute2)});
            family.push_back({childRoute, childFitness});
            
            // Ordenar la familia según la aptitud (de menor a mayor costo).
            sort(family.begin(), family.end(), [](const Individual &a, const Individual &b) {
                return a.fitness < b.fitness;
            });
            
            // Los dos mejores pasan a la siguiente generación.
            newPopulation.push_back(family[0]);
            // Sólo agregar el segundo si aún hay espacio.
            if(newPopulation.size() < populationSize)
                newPopulation.push_back(family[1]);
        }
        
        // Inyección de rutas aleatorias: para cada individuo, con probabilidad pm, se genera una nueva ruta.
        for (auto &indiv : newPopulation) {
            uniform_real_distribution<double> distProb(0.0, 1.0);
            if(distProb(rng) < pm) {
                indiv.route = generateRandomRoute(numCities);
                removeAbruptChanges(indiv.route);
                indiv.fitness = evaluateFitness(indiv.route);
            }
        }
        
        // Actualizar población.
        population = newPopulation;
        
        // (Opcional) Mostrar el mejor individuo de la generación.
        auto bestIt = min_element(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
            return a.fitness < b.fitness;
        });
        cout << "Generación " << gen + 1 << " - Mejor aptitud: " << bestIt->fitness << "\n";
    }
    
    // Mostrar la mejor ruta encontrada.
    auto bestIt = min_element(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
        return a.fitness < b.fitness;
    });
    cout << "\nMejor ruta encontrada:\n";
    for (int cityIndex : bestIt->route) {
        cout << cities[cityIndex].name << " -> ";
    }
    cout << cities[bestIt->route[0]].name << "\n"; // cerrar el circuito (si se desea)
    cout << "Aptitud: " << bestIt->fitness << "\n";
    
    return 0;
}
