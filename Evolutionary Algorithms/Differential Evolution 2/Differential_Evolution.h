#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include <set>
#include <cassert>
#include <sstream>
#include <ctime>

class Solution{
    public:
        Solution () = default;
        Solution (size_t vector_size, double low, double high);
        Solution (std::vector<double> v, double low, double high);
        
        std::string to_string () const;
        double value () const;
        std::vector<double> get_vector() const {return _vector;}
    
    private:
        std::vector<double> _vector;
};


void Differential_Evolution(size_t number_of_solutions,
                            size_t elements_in_vector,
                            double low_bound,
                            double high_bound,
                            double F,
                            int iterations,
                            double crossover_probability);

double random_double (double low, double high);
std::set<int> random_distinct_numbers (int upper_limit, size_t amount_numbers, int except);
std::vector<double> get_donor_vector(Solution const& s1, Solution const& s2, Solution const& s3, double F);
std::vector<double> get_trial_vector(double crossosver_probability, 
                                    std::vector<double> const& original,
                                    std::vector<double> const& donor,
                                    double low,
                                    double high);