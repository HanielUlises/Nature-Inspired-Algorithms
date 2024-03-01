#include <bits/stdc++.h>

void f(double x1, double x2, std::vector<double>& z){
    double goal = (100 * std::pow((x2 - std::pow(x1,2)),2)) + std::pow(x1 - 1, 2);
    z.push_back(goal);
}

double mean (std::vector<double> z){
    double sum = 0;
    for(int i = 0; i < z.size(); i++){
        sum += z[i];
    }
    return sum/z.size();
}

void expected_val (double z, double mean_z, std::vector<double>& exp_z){
    double result = 0;
    result = z / mean_z;
    exp_z.push_back(result);
}

int main (){
    std::vector<std::pair<double, double>> points;

    points.push_back(std::make_pair(-5.5,2.3));
    points.push_back(std::make_pair(1.2,-1.1));
    points.push_back(std::make_pair(6.0,3.0));
    points.push_back(std::make_pair(1.0,0.9));
    points.push_back(std::make_pair(-0.9,0.8));
    points.push_back(std::make_pair(2,4.0));

    // Goal
    std::vector<double> z;
    // Expected value
    std::vector<double> expected_z;

    std::cout << "Points tested" << std::endl;
    for (const auto& point : points) {
        std::cout << "(" << point.first << ", " << point.second << ")" << std::endl;
        f(point.first, point.second, z);
    }

    double mean_z = mean(z);

    for(int i = 0; i < z.size(); i++){
        expected_val(z[i], mean_z, expected_z);
    }

    std::cout<<std::endl<<"Expected values"<<std::endl;

    for(double i: expected_z){
        std::cout<< i<<std::endl;
    }
}