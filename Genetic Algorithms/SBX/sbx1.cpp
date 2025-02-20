#include <bits/stdc++.h>

using std::vector;
using std::pair;
using std::cout;

// Individual, set of variables, both fathers f1 and f2
pair<vector<int>, vector<int>> sbx(vector<double> f1, vector<double> f2){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.00, 1.0);
    double u = dis(gen);
    double beta = 0.0f;
    double beta_c = 0.0f;
    double alpha = 0.0f;
    double nc = 2.0f;

    vector<int> low_bounds (f1.size());
    vector<int> high_bounds (f1.size());

    vector<int> s1 (f1.size());
    vector<int> s2 (f2.size());

    for(int i = 0; i < f1.size() && i < f2.size() ; i++){
        int low_bound = std::min(f1[i], f2[i]);
        int high_bound = std::max(f1[i], f2[i]);

        low_bounds.push_back(low_bound);
        high_bounds.push_back(high_bound);


        alpha = 2 - std::pow(std::abs(beta), -nc + 1);
        beta = 1 + 2/(f1[i] - f2[i]) * std::min((f1[i] - low_bounds[i]), 
                                                (f2[i] - high_bounds[i]));

        beta_c = u <= 1/alpha ? std::pow(u * alpha, 1/(nc+1)) : 
                std::pow(1/(2 - u *alpha),(1/nc)+1) ;
        
        s1[i] = 0.5 * (f1[i] + f2[i] - beta_c * std::abs(f2[i] - f1[i]));
        s2[i] = 0.5 * (f1[i] + f2[i] + beta_c * std::abs(f2[i] - f1[i]));
    }
    pair<vector<int>, vector<int>> offspring = std::make_pair(s1, s2);
    return offspring;
}

int main (){
    vector<double> f1 = {2.3, 4.5};
    vector<double> f2 = {1.4, -0.2};
    pair<vector<int>, vector<int>> off = sbx(f1, f2);
}