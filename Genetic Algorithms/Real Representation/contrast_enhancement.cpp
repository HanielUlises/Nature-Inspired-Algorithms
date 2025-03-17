#include <opencv2/opencv.hpp>
#include <iostream>
#include "genetic_algo_img.h"
#include "image_utils.h"

int main() {
    std::string path, chosen_image;
    cv::Mat gray_img;

    std::cout << "Select image:\n";
    std::cout << "1. Original_Medica5R\n";
    std::cout << "2. Original_Medica6R\n";
    std::cout << "Enter choice (1 or 2): ";

    int choice;
    std::cin >> choice;

    if (choice == 1) {
        path = "images/Original_Medica5R.png";
        chosen_image = "First";
        gray_img = cv::imread(path, cv::IMREAD_GRAYSCALE);
    } else if (choice == 2) {
        path = "images/Original_Medica6R.png";
        chosen_image = "Second";
        gray_img = cv::imread(path, cv::IMREAD_GRAYSCALE);
    } else {
        std::cout << "Invalid image choice. Exiting.\n";
        return -1;
    }

    if (gray_img.empty()) {
        std::cout << "Image file not found: " << path << std::endl;
        return -1;
    }
    gray_img.convertTo(gray_img, CV_32F, 1.0 / 255.0);

    std::cout << "Select fitness metric:\n";
    std::cout << "1. Shannon Entropy\n";
    std::cout << "2. Standard Deviation\n";
    std::cout << "Enter choice (1 or 2): ";

    std::cin >> choice;

    GeneticAlgorithmReal::FitnessMetric metric;
    std::string metric_name;
    if (choice == 1) {
        metric = GeneticAlgorithmReal::FitnessMetric::ENTROPY;
        metric_name = "Entropy";
    } else if (choice == 2) {
        metric = GeneticAlgorithmReal::FitnessMetric::STDDEV;
        metric_name = "Standard Deviation";
    } else {
        std::cout << "Invalid metric choice. Defaulting to Shannon Entropy.\n";
        metric = GeneticAlgorithmReal::FitnessMetric::ENTROPY;
        metric_name = "Entropy";
    }

    int pop_size = 20;
    int num_genes = 2;
    int num_generations = 50;
    double crossover_prob = 0.8;
    double mutation_prob = 0.1;
    std::vector<double> lower_bound = {0.0, 0.0};
    std::vector<double> upper_bound = {20.0, 1.0};

    GeneticAlgorithmReal ga(pop_size, num_genes, num_generations, crossover_prob, mutation_prob,
                            lower_bound, upper_bound, gray_img, metric);
    ga();

    auto [best_solution, best_fitness] = ga.getBestSolution();
    double alpha = best_solution[0];
    double delta = best_solution[1];
    std::cout << "Best alpha: " << alpha << ", Best delta: " << delta
              << ", " << metric_name << ": " << -best_fitness << std::endl;

    cv::Mat new_image = sigmoid_transform(gray_img, alpha, delta);
    cv::normalize(new_image, new_image, 0, 1, cv::NORM_MINMAX);
    new_image.convertTo(new_image, CV_8U, 255.0);

    std::string final = "transformed_image_" + chosen_image + "_" + metric_name + ".png";
    cv::imwrite(final, new_image);

    std::cout << "Transformed image saved as: " << final << std::endl;

    return 0;
}