#include "tsp_solver.h"
#include <iostream>
#include <limits>
#include <iomanip>
#include <string>

#define INFINITY std::numeric_limits<double>::infinity()

// #define TIME_WINDOWS
#define NO_TIME

int main() {
    std::vector<std::vector<double>> distances = {
        {0,    61.82, 18.54, 37.52, 54.08,  1.88, 59.98, 32.82, 69.42, 36.76, 60.26},
        {61.82, 0,    50.84, 33.62,  7.50, 59.88,  2.76, 28.84,  7.78, 28.14,  5.80},
        {18.54,50.84,  0,    26.74, 43.38, 18.60, 49.28, 22.00, 58.70, 23.36, 49.30},
        {37.52,33.62, 26.74,  0,    26.16, 35.56, 32.06,  4.80, 41.50,  3.26, 32.08},
        {54.08, 7.50, 43.38, 26.16,  0,    52.06,  7.32, 21.38, 15.34, 20.68,  5.92},
        { 1.88,59.88, 18.60, 35.56, 52.06,  0,    57.96, 30.86, 67.38, 34.80, 58.30},
        {59.98, 2.76, 49.28, 32.06,  7.32, 57.96,  0,    27.28, 10.62, 26.58,  6.76},
        {32.82,28.84, 22.00,  4.80, 21.38, 30.86, 27.28,  0,    36.72,  4.02, 27.30},
        {69.42, 7.78, 58.70, 41.50, 15.34, 67.38, 10.62, 36.72,  0,    36.02, 12.14},
        {36.76,28.14, 23.36,  3.26, 20.68, 34.80, 26.58,  4.02, 36.02,  0,    26.60},
        {60.26, 5.80, 49.30, 32.08,  5.92, 58.30,  6.76, 27.30, 12.14, 26.60,  0   }
    };

#ifdef TIME_WINDOWS
    std::vector<City> cities = {
        {"New York",    0, INFINITY},
        {"Los Angeles", 50, 90},
        {"Chicago",    15, 25},
        {"Houston",    30, 55},
        {"Phoenix",    15, 75},
        {"Philadelphia",5, 35},
        {"San Diego", 150, 200},
        {"Dallas",     25, 50},
        {"San Francisco", 65, 100},
        {"Austin",    120, 150},
        {"Las Vegas",  30, 85}
    };
#endif

#ifdef NO_TIME
    std::vector<City> cities = {
        {"New York",    -INFINITY, INFINITY},
        {"Los Angeles", -INFINITY, INFINITY},
        {"Chicago",    -INFINITY, INFINITY},
        {"Houston",    -INFINITY, INFINITY},
        {"Phoenix",    -INFINITY, INFINITY},
        {"Philadelphia",-INFINITY, INFINITY},
        {"San Diego", -INFINITY, INFINITY},
        {"Dallas",     -INFINITY, INFINITY},
        {"San Francisco", -INFINITY, INFINITY},
        {"Austin",    -INFINITY, INFINITY},
        {"Las Vegas", -INFINITY, INFINITY}
    };
#endif

    if (distances.size() != cities.size() || distances.empty()) {
        std::cerr << "Error: La matriz de distancias y el vector de ciudades no coinciden en tamaño o están vacíos.\n";
        return 1;
    }
    for (const auto& row : distances) {
        if (row.size() != distances.size()) {
            std::cerr << "Error: La matriz de distancias no es cuadrada.\n";
            return 1;
        }
    }

    TspSolver solver(distances, cities,
                     /*pop_size=*/30,
                     /*generations=*/100,
                     /*cross_rate=*/0.9,
                     /*mut_rate=*/0.15,
                     /*rand_insert_rate=*/0.05,
                     /*penalty=*/10.0);

    solver.solve();

    std::cout << "\n==== Resultado Final ====\n";
    std::cout << "Distancia total de la mejor ruta: "
              << std::fixed << std::setprecision(2) << solver.get_best_distance() << "\n\n";
    std::cout << "Ruta óptima:\n";

    const auto& best_tour = solver.get_best_tour();
    double distancia_acum = 0.0;

    std::cout << std::left << std::setw(4) << "Paso"
              << std::setw(20) << "Ciudad"
              << std::setw(12) << "Distancia (d)"
              << "\n";
    std::cout << std::string(36, '-') << "\n";

    for (size_t i = 0; i < best_tour.size(); ++i) {
        int curr = best_tour[i];
        std::cout << std::left << std::setw(4) << i + 1
                  << std::setw(20) << cities[curr].name;

        if (i > 0) {
            int prev = best_tour[i - 1];
            double d = distances[prev][curr];
            distancia_acum += d;
            std::cout << std::setw(12) << std::fixed << std::setprecision(2) << distancia_acum;
        } else {
            std::cout << std::setw(12) << "0.00";
        }
        std::cout << "\n";
    }

    if (!best_tour.empty()) {
        int last = best_tour.back();
        int origin = best_tour.front();
        double d_back = distances[last][origin];
        distancia_acum += d_back;

        std::cout << std::string(36, '-') << "\n";
        std::cout << std::left << std::setw(4) << best_tour.size() + 1
                  << std::setw(20) << cities[origin].name + " (regreso)"
                  << std::setw(12) << std::fixed << std::setprecision(2) << distancia_acum
                  << "\n";
    }

    return 0;
}