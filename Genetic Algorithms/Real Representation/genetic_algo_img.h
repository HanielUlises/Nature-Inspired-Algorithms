#ifndef GENETIC_ALGO_IMG_H
#define GENETIC_ALGO_IMG_H

#include <vector>
#include <random>
#include <utility>
#include <fstream>
#include <opencv2/opencv.hpp>
#include "image_utils.h"
#include "PM/polynomial_mut.h"
#include "SBX/sbx.h"

class GeneticAlgorithmReal {
public:
    enum class FitnessMetric { ENTROPY, STDDEV };

private:
    int population_size;
    int num_generations;
    int num_genes;
    double crossover_prob;
    double mutation_prob;
    std::vector<double> lower_bound;
    std::vector<double> upper_bound;
    cv::Mat gray_img;
    std::vector<std::vector<double>> population;
    std::vector<double> fitness_values;
    std::ofstream output_file;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    FitnessMetric metric;

    void initialize_population();
    void evaluate_fitness();
    int tournament_selection();
    void crossover();
    void mutate();
    void select_best_individual();
    void apply_elitism();

public:
    GeneticAlgorithmReal(int pop_size, int num_genes, int num_generations,
                         double crossover_prob, double mutation_prob,
                         std::vector<double> lower_bound,
                         std::vector<double> upper_bound,
                         const cv::Mat& gray_img,
                         FitnessMetric metric);
    ~GeneticAlgorithmReal();
    void operator()();
    std::pair<std::vector<double>, double> getBestSolution();
};

#endif // GENETIC_ALGO_IMG_H