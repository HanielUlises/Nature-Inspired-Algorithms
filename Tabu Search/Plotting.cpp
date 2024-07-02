#include "matplotlibcpp.h"
#include "tabu_search.h"

namespace plt = matplotlibcpp;

void TabuSearch::plotConvergenceGraph() {
    

    std::string title = "Convergence Graph ";
    title.append("Tabu Search");

    std::vector<int> generations(evaluation_history.size());
    std::iota(generations.begin(), generations.end(), 0);

    plt::plot(generations, evaluation_history, {{"label", "evaluations"}});

    plt::xlabel("Generation");
    plt::ylabel("Fitness");
    plt::legend();

    plt::title(title);

    plt::show();
}
