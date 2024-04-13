#include "matplotlibcpp.h"
#include <string>
#include <fstream>
#include "ploteado.h"

void plotConvergenceGraph(std::vector<double> average,std::vector<double> best, std::vector<double> worst) {
    namespace plt = matplotlibcpp;

    std::string title = "Convergence Graph ";

    std::vector<int> generations(best.size());
    std::iota(generations.begin(), generations.end(), 0);

    plt::plot(generations, best, {{"label", "best"}});
    plt::plot(generations, average, {{"label", "worst"}});
    plt::plot(generations, worst, {{"label", "average"}});

    plt::xlabel("Generation");
    plt::ylabel("Fitness");
    plt::legend();

    plt::title(title);

    plt::show();
}

