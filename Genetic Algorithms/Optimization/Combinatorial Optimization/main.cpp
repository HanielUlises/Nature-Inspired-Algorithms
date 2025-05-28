#include <iostream>
#include <vector>
#include <limits>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include "tsp_solver.h"

void print_separator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "  " << title << std::endl;
    std::cout << std::string(60, '=') << std::endl;
}

void print_route_details(const std::vector<int>& route, const std::vector<City>& cities, 
                        const std::vector<std::vector<double>>& distances, double total_time) {
    if (route.empty()) {
        std::cout << "Ruta vacía" << std::endl;
        return;
    }
    
    std::cout << "Ruta: ";
    for (size_t i = 0; i < route.size(); ++i) {
        std::cout << cities[route[i]].name;
        if (i < route.size() - 1) std::cout << " -> ";
    }
    std::cout << " -> " << cities[route[0]].name << std::endl;
    
    std::cout << "Tiempo total: " << std::fixed << std::setprecision(2) << total_time << std::endl;
    
    double current_time = 0.0;
    std::cout << "\nDetalles de la ruta:" << std::endl;
    std::cout << "Ciudad de inicio: " << cities[route[0]].name << " (tiempo: 0.00)" << std::endl;
    
    for (size_t i = 0; i < route.size() - 1; ++i) {
        int from = route[i];
        int to = route[i + 1];
        double travel_time = distances[from][to];
        current_time += travel_time;
        
        double arrival_time = current_time;
        if (current_time < cities[to].window_start) {
            current_time = cities[to].window_start;
        }
        
        std::cout << cities[from].name << " -> " << cities[to].name 
                  << " (viaje: " << std::fixed << std::setprecision(2) << travel_time
                  << ", llegada: " << arrival_time;
        
        if (arrival_time != current_time) {
            std::cout << ", espera hasta: " << current_time;
        }
        
        std::cout << ", ventana: [" << cities[to].window_start << ", " << cities[to].window_end << "]";
        
        if (arrival_time > cities[to].window_end) {
            std::cout << " *** VIOLACIÓN ***";
        }
        
        std::cout << ")" << std::endl;
    }
    
    double return_time = distances[route.back()][route[0]];
    current_time += return_time;
    std::cout << cities[route.back()].name << " -> " << cities[route[0]].name 
              << " (regreso: " << std::fixed << std::setprecision(2) << return_time
              << ", tiempo final: " << current_time << ")" << std::endl;
}

void run_experiment(const std::string& experiment_name, 
                   const std::vector<City>& cities,
                   const std::vector<std::vector<double>>& distances,
                   int num_runs = 5, int generations = 100) {
    
    print_separator(experiment_name);
    
    std::vector<double> times;
    std::vector<std::vector<int>> routes;
    double best_overall_time = std::numeric_limits<double>::infinity();
    std::vector<int> best_overall_route;
    
    for (int run = 1; run <= num_runs; ++run) {
        std::cout << "\n--- Ejecución " << run << " ---" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        TspSolver solver(cities, distances);
        solver.run(100, generations, 0.1, 0.05);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::vector<int> route = solver.get_best_route();
        double time = solver.get_best_route_time();
        
        times.push_back(time);
        routes.push_back(route);
        
        std::cout << "Tiempo de ejecución: " << duration.count() << " ms" << std::endl;
        
        if (std::isinf(time)) {
            std::cout << "No se encontró ruta válida (inf)" << std::endl;
        } else {
            std::cout << "Mejor tiempo encontrado: " << std::fixed << std::setprecision(2) << time << std::endl;
            if (time < best_overall_time) {
                best_overall_time = time;
                best_overall_route = route;
            }
        }
    }
    
    std::cout << "\n" << std::string(40, '-') << std::endl;
    std::cout << "RESUMEN ESTADÍSTICO" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    std::vector<double> valid_times;
    for (double t : times) {
        if (!std::isinf(t)) {
            valid_times.push_back(t);
        }
    }
    
    std::cout << "Ejecuciones válidas: " << valid_times.size() << "/" << num_runs << std::endl;
    
    if (!valid_times.empty()) {
        double sum = 0.0;
        double min_time = *std::min_element(valid_times.begin(), valid_times.end());
        double max_time = *std::max_element(valid_times.begin(), valid_times.end());
        
        for (double t : valid_times) {
            sum += t;
        }
        double avg_time = sum / valid_times.size();
        
        std::cout << "Mejor tiempo: " << std::fixed << std::setprecision(2) << min_time << std::endl;
        std::cout << "Peor tiempo: " << std::fixed << std::setprecision(2) << max_time << std::endl;
        std::cout << "Tiempo promedio: " << std::fixed << std::setprecision(2) << avg_time << std::endl;
        
        if (!best_overall_route.empty()) {
            std::cout << "\n--- MEJOR RUTA ENCONTRADA ---" << std::endl;
            print_route_details(best_overall_route, cities, distances, best_overall_time);
        }
    } else {
        std::cout << "No se encontraron rutas válidas en ninguna ejecución." << std::endl;
    }
    
    std::cout << "\nResultados de todas las ejecuciones:" << std::endl;
    for (int i = 0; i < num_runs; ++i) {
        std::cout << "Ejecución " << (i+1) << ": ";
        if (std::isinf(times[i])) {
            std::cout << "inf (sin solución válida)" << std::endl;
        } else {
            std::cout << std::fixed << std::setprecision(2) << times[i] << std::endl;
        }
    }
}

int main() {
    std::vector<std::vector<double>> distances = {
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

    std::vector<City> cities_with_windows = {
        {"New York", 0, std::numeric_limits<double>::infinity()},  
        {"Los Angeles", 50, 90},
        {"Chicago", 15, 25},
        {"Houston", 30, 55},
        {"Phoenix", 15, 75},
        {"Philadelphia", 5, 35},
        {"San Diego", 150, 200},
        {"Dallas", 25, 50},
        {"San Francisco", 65, 100},
        {"Austin", 120, 150},
        {"Las Vegas", 30, 85}
    };

    std::vector<City> cities_no_windows = {
        {"New York", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Los Angeles", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Chicago", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Houston", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Phoenix", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Philadelphia", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"San Diego", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Dallas", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"San Francisco", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Austin", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {"Las Vegas", -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}
    };

    std::cout << "ALGORITMO GENÉTICO PARA TSP CON VENTANAS DE TIEMPO" << std::endl;
    std::cout << "Parámetros: Población=100, Generaciones=100, Mutación=0.1, Inyección=0.05" << std::endl;

    run_experiment("EXPERIMENTO 1: CON VENTANAS DE TIEMPO", cities_with_windows, distances, 5, 100);
    run_experiment("EXPERIMENTO 2: SIN VENTANAS DE TIEMPO", cities_no_windows, distances, 5, 100);

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "  EXPERIMENTOS COMPLETADOS" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    return 0;
}