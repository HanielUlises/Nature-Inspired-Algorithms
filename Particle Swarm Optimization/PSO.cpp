#include "PSO.h"
#include "LinealRegression.h"
#include <cmath>
#include <limits>
#include <random>
#include <iostream>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

Particle::Particle(int dimensions) : best_score(std::numeric_limits<double>::infinity()) {
    position.resize(dimensions);
    velocity.resize(dimensions);
    best_position.resize(dimensions);

    for (int i = 0; i < dimensions; ++i) {
        position[i] = dis(gen) * 10 - 5; // Initialize to random values within a range
        velocity[i] = dis(gen) * 2 - 1; // Random velocities
    }
    best_position = position;
}

void Particle::updateVelocity(const std::vector<double>& global_best_position, double omega, double phi_p, double phi_g) {
    for (size_t i = 0; i < position.size(); ++i) {
        double r_p = dis(gen);
        double r_g = dis(gen);
        velocity[i] = omega * velocity[i] +
                      phi_p * r_p * (best_position[i] - position[i]) +
                      phi_g * r_g * (global_best_position[i] - position[i]);
    }
}

void Particle::updatePosition() {
    for (size_t i = 0; i < position.size(); ++i) {
        position[i] += velocity[i];
    }
}

PSO::PSO(int swarm_size, int dimensions, std::function<double(const std::vector<double>&)> obj_function, int max_iterations)
  : global_best_score(std::numeric_limits<double>::infinity()),
    objective_function(obj_function) { 
    particles = std::vector<Particle>(swarm_size, Particle(dimensions));
    global_best_position.resize(dimensions);
    history_global_best_score.resize(max_iterations);
}

void PSO::optimize(int max_iterations, double omega, double phi_p, double phi_g) {
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        for (auto& p : particles) {
            p.updateVelocity(global_best_position, omega, phi_p, phi_g);
            p.updatePosition();
            double score = objective_function(p.position);
            if (score < p.best_score) {
                p.best_score = score;
                p.best_position = p.position;
            }
            if (score < global_best_score) {
                global_best_score = score;
                global_best_position = p.best_position;
            }
        }
        history_global_best_score.push_back(global_best_score);
    }
    
}

void PSO::printResults() const {
    std::cout << "Best position: ";
    for (auto x : global_best_position) {
        std::cout << x << " ";
    }
    std::cout << "\nBest score: " << global_best_score << std::endl;
}

PSO::PSO(int swarm_size, int dimensions)
  : global_best_score(std::numeric_limits<double>::infinity()){ 
    particles = std::vector<Particle>(swarm_size, Particle(dimensions));
    global_best_position.resize(dimensions);
}

void PSO::minimizeError(int max_iterations, double omega, double phi_p, double phi_g) {
    std::vector<std::vector<double>> DataSet = readCSV();
    plottingSLR(DataSet);
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        for (auto& p : particles) {
            p.updateVelocity(global_best_position, omega, phi_p, phi_g);
            p.updatePosition();
            double score = evaluation(p.position,DataSet);
            if (score < p.best_score) {
                p.best_score = score;
                p.best_position = p.position;
            }
            if (score < global_best_score) {
                global_best_score = score;
                global_best_position = p.best_position;
            }
        }
    }
    plottingSLR_withSolution(DataSet, global_best_position);
    

}

