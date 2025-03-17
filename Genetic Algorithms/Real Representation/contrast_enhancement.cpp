#include <opencv2/opencv.hpp>
#include <iostream>
#include "genetic_algo_img.h"
#include "image_utils.h"

int main() {
    std::string path = "images/Original_Medica6R.png";
    cv::Mat gray_img = cv::imread(path, cv::IMREAD_GRAYSCALE);
    if (gray_img.empty()) {
        std::cout << "Image file not found" << std::endl;
        return -1;
    }
    gray_img.convertTo(gray_img, CV_32F, 1.0 / 255.0);


    int pop_size = 20;
    int num_genes = 2; 
    int num_generations = 50;
    double crossover_prob = 0.8;
    double mutation_prob = 0.1;
    std::vector<double> lower_bound = {0.0, 0.0};    
    std::vector<double> upper_bound = {20.0, 1.0};   

    GeneticAlgorithmReal ga(pop_size, num_genes, num_generations, crossover_prob, mutation_prob,
                            lower_bound, upper_bound, gray_img);
    ga(); 

    auto [best_solution, best_fitness] = ga.getBestSolution();
    double alpha = best_solution[0];
    double delta = best_solution[1];
    std::cout << "Best alpha: " << alpha << ", Best delta: " << delta 
              << ", Entropy: " << -best_fitness << std::endl;

    cv::Mat new_image = sigmoid_transform(gray_img, alpha, delta);
    cv::normalize(new_image, new_image, 0, 1, cv::NORM_MINMAX);
    new_image.convertTo(new_image, CV_8U, 255.0);
    cv::imwrite("transformed_image_2.png", new_image);

    return 0;
}