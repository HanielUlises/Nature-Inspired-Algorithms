#include "Plotting.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

void plotting(std::vector<std::vector<double>> history){
    std::vector<double> bestscore=history[0];
    std::vector<double> worstscore=history[0];
    for (const auto &score : history){
        if (score[score.size()-1]<=bestscore[bestscore.size()-1]) bestscore=score;
    }
    for (const auto &score : history){
        if (score[score.size()-1]>=worstscore[worstscore.size()-1]) worstscore=score;
    }
    std::string title = "Convergence Graph ";
    //title.append(function);

    std::vector<int> generations(bestscore.size());
    std::iota(generations.begin(), generations.end(), 0);

    plt::plot(generations, bestscore, {{"label", "best"}});
    plt::plot(generations, worstscore, {{"label", "worst"}});

    plt::xlabel("Generation");
    plt::ylabel("Best global score");
    plt::legend();

    plt::title(title);

    plt::show();
    
    
}