#include <numeric>
#include <vector>
#include <utility>

#define M_PI 3.14159265358979323846
#define M_E 2.71828182845904523536

double rosenbrockFunction(const std::vector<double>& individual);
double ackleyFunction(const std::vector<double>& individual);
double GriewankFunction(const std::vector<double>& individual);
double RastriginFunction(const std::vector<double>& individual);
double SimpleLinearRegression(const std::vector<double>& individual);
double Langermann(const std::vector<double> &individual);
double dropWave(const std::vector<double> &individual);