#include <bits/stdc++.h>

using std::vector;
using std::pair;
using std::cout;

// SBX Crossover
pair<vector<double>, vector<double>> sbx(vector<double> f1, vector<double> f2){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double nc = 2.0f; 

    vector<double> low_bounds(f1.size());
    vector<double> high_bounds(f1.size());

    vector<double> s1(f1.size());
    vector<double> s2(f2.size());

    for(size_t i = 0; i < f1.size(); i++){
        double low_bound = std::min(f1[i], f2[i]);
        double high_bound = std::max(f1[i], f2[i]);

        low_bounds[i] = low_bound;
        high_bounds[i] = high_bound;

        if (std::abs(f1[i] - f2[i]) < 1e-14) {
            s1[i] = f1[i];
            s2[i] = f2[i];
            continue;
        }

        double beta = 1 + (2.0 * std::min(f1[i] - low_bound, high_bound - f2[i]) / 
                           std::abs(f1[i] - f2[i]));

        double alpha = 2 - std::pow(beta, -(nc + 1));

        double u = dis(gen);
        double beta_c = (u <= 1.0 / alpha) ? 
                        std::pow((u * alpha), (1.0 / (nc + 1))) : 
                        std::pow((1.0 / (2.0 - u * alpha)), (1.0 / (nc + 1)));

        s1[i] = 0.5 * (f1[i] + f2[i] - beta_c * std::abs(f2[i] - f1[i]));
        s2[i] = 0.5 * (f1[i] + f2[i] + beta_c * std::abs(f2[i] - f1[i]));
    }

    return {s1, s2};
}

int main() {
    vector<double> f1 = {2.3, 4.5};
    vector<double> f2 = {1.4, -0.2};

    pair<vector<double>, vector<double>> off = sbx(f1, f2);

    cout << "Offspring 1: ";
    for(double v : off.first) cout << v << " ";
    cout << "\nOffspring 2: ";
    for(double v : off.second) cout << v << " ";
    cout << "\n";

    return 0;
}
