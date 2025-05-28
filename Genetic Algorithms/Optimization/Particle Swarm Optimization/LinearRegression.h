#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H

#include <vector>

std::vector<std::vector<double>> readCSV();
void plottingSLR(std::vector<std::vector<double>> DataSet);
void plotting(const std::vector<std::vector<double>>& history) ;
double evaluation(std::vector<double> ab, std::vector<std::vector<double>> DataSet);
void plottingSLR_withSolution(std::vector<std::vector<double>> DataSet, std::vector<double> particle);

#endif