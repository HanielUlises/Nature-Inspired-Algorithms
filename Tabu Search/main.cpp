#include "tabu_search.h"
#include <fstream>
#include <sstream>

std::vector<Object> readObjectsFromFile(const std::string& filename, double& container_capacity) {
    std::ifstream file(filename);
    std::string line;
    std::vector<Object> objects;
    int id = 1;

    if (file.is_open()) {
        std::getline(file, line);
        std::istringstream ss(line);
        ss >> container_capacity;

        while (std::getline(file, line)) {
            double weight;
            std::istringstream obj_ss(line);
            obj_ss >> weight;
            objects.push_back({id++, weight});
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
    return objects;
}

int main() {
    double container_capacity;
    std::vector<Object> objects = readObjectsFromFile("input.txt", container_capacity);

    int max_iterations = 1000;
    int tabu_list_size = 50;

    TabuSearch tabu_search(objects, container_capacity, max_iterations, tabu_list_size);
    tabu_search.solve();
    tabu_search.printSolution();

    return 0;
}
