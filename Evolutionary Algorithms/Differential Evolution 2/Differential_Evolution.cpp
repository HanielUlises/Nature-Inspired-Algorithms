#include "Differential_Evolution.h"

double random_double (double low, double high){
    double d = (double) rand() / RAND_MAX;
    return low + d * (high- low);
}

std::set<int> random_distinct_numbers (int upper_limit, size_t amount_numbers, int except){
    std::set <int> numbers;
    while(numbers.size() < amount_numbers){
        int r = rand() % upper_limit;
        if(r != except){
            numbers.insert(r);
        }
    }
    return numbers;
}

Solution::Solution(size_t vector_size, double low, double high){
    _vector.resize(vector_size);
    std::generate_n(_vector.begin(), vector_size, [low, high]() {return random_double(low, high);});
}

Solution::Solution(std::vector<double> v, double low, double high) : _vector{v} {}

std::string Solution::to_string() const{
    std::ostringstream stream;
    stream << '[]';

    for(double d: _vector){
        stream << d << ' ';
    }

    stream << "] v = " << value() << std::endl;
    return stream.str();
}

double Solution::value() const{
    return 0;
}

std::vector<double> get_donor_vector(Solution const& s1, Solution const& s2, Solution const& s3, double F){
    std::vector<double> dv;
    size_t size = s2.get_vector().size();

    for(size_t i = 0; i < size; i++){
        dv.push_back(s1.get_vector()[i] + F * (s2.get_vector()[i] - s3.get_vector()[i]));
    }
    return dv;
}

std::vector<double> get_trial_vector(double crossosver_probability, 
                                    std::vector<double> const& original,
                                    std::vector<double> const& donor,
                                    double low,
                                    double high){
    std::vector<double> trial;
    int stable = rand() % original.size();

    for(size_t i = 0; i < original.size(); i++){
        if(i == stable){
            trial.push_back(donor[i]);
        }else{
            double r = rand() / RAND_MAX;

            if(r < crossosver_probability){
                trial.push_back(donor[i]);
            }else{
                trial.push_back(original[i]);
            }
        }
    }
    return trial;
}

void Differential_Evolution(size_t number_of_solutions,
                            size_t elements_in_vector,
                            double low_bound,
                            double high_bound,
                            double F,
                            int iterations,
                            double crossover_probability){

   assert(number_of_solutions >= 4);

   std::vector<Solution> Solutions;

   Solutions.resize(number_of_solutions);

   std::generate_n(Solutions.begin(), number_of_solutions, [elements_in_vector, low_bound, high_bound](){
    Solution s(elements_in_vector, low_bound, high_bound);
        return s;
   });

   for(int it = 0; it < iterations; it++){
    std::vector<Solution> trial_solutions;
    for(size_t i = 0; i < number_of_solutions; i++){
        std::set<int> indexes = random_distinct_numbers(number_of_solutions, 3, i);
        std::vector<Solution> donor_components;

        for(int index : indexes){
            donor_components.push_back(Solutions[index]);
        }

        std::vector<double> donor_vector = get_donor_vector(donor_components[0],donor_components[1], donor_components[2], F);
        std::vector<double> trial_vector = get_trial_vector(crossover_probability, Solutions[i].get_vector(), donor_vector, low_bound, high_bound);

        Solution trial_solution (trial_vector, low_bound, high_bound);
        trial_solutions.push_back(trial_solution);
    }
    for(size_t i = 0; i < number_of_solutions; i++){
        if(trial_solutions[i].value() < Solutions[i].value()){
            Solutions[i] = trial_solutions[i];
        }
        std::cout<< Solutions[i].to_string() << std::endl;
    }
   }                     
}