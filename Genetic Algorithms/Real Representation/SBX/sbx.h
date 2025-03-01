#ifndef SBX_H
#define SBX_H

#include <vector>
#include <utility>
#include <cmath>
#include <random>

std::pair<std::vector<double>, std::vector<double>> sbx(std::vector<double> f1, std::vector<double> f2, double r, 
    const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, double nc = 2.0);
    
#endif // SBX_H
