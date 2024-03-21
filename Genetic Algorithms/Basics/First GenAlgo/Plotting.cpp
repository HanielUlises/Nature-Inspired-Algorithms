// Had to make a plotting file in order to solve reference problems of matplotlibcpp
#include "matplotlibcpp.h"
#include "GeneticAlgorithm.h"

namespace plt = matplotlibcpp;

void GeneticAlgorithm::plotConvergenceGraph(std::string function) {
    

    std::string title = "Convergence Graph ";
    title.append(function);

    std::vector<int> generations(bestFitnessHistory.size());
    std::iota(generations.begin(), generations.end(), 0);

    plt::plot(generations, bestFitnessHistory, {{"label", "best"}});
    plt::plot(generations, worstFitnessHistory, {{"label", "worst"}});
    plt::plot(generations, averageFitnessHistory, {{"label", "average"}});

    plt::xlabel("Generation");
    plt::ylabel("Fitness");
    plt::legend();

    plt::title(title);

    plt::show();
}
