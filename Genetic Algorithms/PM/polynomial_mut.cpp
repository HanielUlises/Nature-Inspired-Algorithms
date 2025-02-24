#include <bits/stdc++.h>

using std::vector;
using std::pair;
using std::cout;

vector<double> polynomial_mutation(vector<double>& x, double r, int nm) {
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<> dis(0.0, 1.0);

    vector<double> mutated_x(x.size());
    
    
    for (int i = 0; i < x.size(); i++) {
        double upper_bound = x[i] > 0 ? floor(x[i]) : ceil(x[i]);
        double lower_bound = x[i] < 0 ? floor(x[i]) : ceil(x[i]);

        double delta = std::min(upper_bound - x[i], x[i] - lower_bound) / (upper_bound - lower_bound);
        double deltaq = 0.0;

        if (r <= 0.5) {
            deltaq = x[i] + std::pow(2 * delta + (1 - 2 * delta) * (1 - r), nm + 1);
        } else {
            deltaq = x[i] - std::pow(2 * (1 - delta) + 2 * (delta - 0.5) * (1 - r), nm + 1);
        }

        mutated_x[i] = deltaq * (upper_bound - lower_bound);
    }
    
    return mutated_x;
}

int main() {
    vector<double> father = {1.126, -0.588};  
    double r = 0.4;
    int nm = 20;  

    vector<double> mutated_offspring = polynomial_mutation(father, r, nm);
    
    cout << "Mutated offspring: ";
    for (const auto& val : mutated_offspring) {
        cout << val << " ";
    }
    cout << std::endl;

    return 0;
}
