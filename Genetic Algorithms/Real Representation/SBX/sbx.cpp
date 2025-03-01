#include "sbx.h"

// SBX Crossover
std::pair<std::vector<double>, std::vector<double>> sbx(std::vector<double> f1, std::vector<double> f2, double r, 
    const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, double nc) {
    std::vector<double> s1(f1.size());
    std::vector<double> s2(f2.size());

    for(size_t i = 0; i < f1.size(); i++) {
        double low_bound_value = f1[i];
        double high_bound_value = f2[i];


        if (std::abs(f1[i] - f2[i]) < 1e-14) {
            s1[i] = f1[i];
            s2[i] = f2[i];
            continue;
        }

        double beta = 1 + (2.0 * std::min(f1[i] - low_bound_value, high_bound_value - f2[i]) / 
                           std::abs(f1[i] - f2[i]));

        double alpha = 2 - std::pow(std::abs(beta), -(2.0 + 1.0));

        double u = r;
        double beta_c = (u <= 1.0 / alpha) ? 
                        std::pow((u * alpha), (1.0 / (2 + 1))) : 
                        std::pow((1.0 / (2.0 - u * alpha)), (1.0 / (2 + 1)));

        s1[i] = 0.5 * (f1[i] + f2[i] - beta_c * std::abs(f2[i] - f1[i]));
        s2[i] = 0.5 * (f1[i] + f2[i] + beta_c * std::abs(f2[i] - f1[i]));
    }

    return {s1, s2};
}