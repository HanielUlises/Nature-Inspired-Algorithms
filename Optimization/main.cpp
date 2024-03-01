#include "function_test_1.h"

int main (){
    // Parte uno
    std::vector<std::pair<double, double>> points_p1;

    points_p1.push_back(std::make_pair(-5.5,2.3));
    points_p1.push_back(std::make_pair(1.2,-1.1));
    points_p1.push_back(std::make_pair(6.0,3.0));
    points_p1.push_back(std::make_pair(1.0,0.9));
    points_p1.push_back(std::make_pair(-0.9,0.8));
    points_p1.push_back(std::make_pair(2,4.0));

    // Goal
    std::vector<double> z;
    // Expected value
    std::vector<double> expected_z;

    std::cout << "Points 1 tested" << std::endl;
    for (const auto& point : points_p1) {
        std::cout << "(" << point.first << ", " << point.second << ")" << std::endl;
        f1(point.first, point.second, z);
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