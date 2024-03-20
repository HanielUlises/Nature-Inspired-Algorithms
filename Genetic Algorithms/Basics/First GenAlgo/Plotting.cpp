// Plotting.cpp
#include "matplotlibcpp.h"
#include "GeneticAlgorithm.h"

void GeneticAlgorithm::plotConvergenceGraph() {
    namespace plt = matplotlibcpp;

    std::vector<int> generations(bestFitnessHistory.size());
    std::iota(generations.begin(), generations.end(), 0);

    plt::plot(generations, bestFitnessHistory, {{"label", "best"}});
    plt::plot(generations, worstFitnessHistory, {{"label", "worst"}});
    plt::plot(generations, averageFitnessHistory, {{"label", "average"}});

    plt::xlabel("Generation");
    plt::ylabel("Fitness");
    plt::legend();
    plt::title("Convergence Graph");

    // Show the plot
    plt::show();
}
