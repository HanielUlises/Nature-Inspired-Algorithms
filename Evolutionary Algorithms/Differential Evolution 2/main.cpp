#include "Differential_Evolution.h"

int main (){
    srand(time(NULL));

    const size_t k_number_of_solutions = 5;
    const size_t k_elements_in_vector = 3;
    const size_t k_number_iterations = 100;

    const double k_low_bound = -10.0f;
    const double k_high_bound = 10.0f;
    const double F = 0.98f;
    const double k_crossosver_probability = 0.8f;

    Differential_Evolution(k_number_of_solutions,
                            k_elements_in_vector,
                            k_low_bound,
                            k_high_bound,
                            F,
                            k_number_iterations, 
                            k_crossosver_probability);

    return 0;
}