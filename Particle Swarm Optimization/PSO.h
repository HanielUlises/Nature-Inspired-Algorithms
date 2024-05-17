#ifndef PSO_H
#define PSO_H

#include <vector>
#include <functional>

class Particle {
public:
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> best_position;
    double best_score;

    Particle(int dimensions);
    void updateVelocity(const std::vector<double>& global_best_position, double omega, double phi_p, double phi_g);
    void updatePosition();
};

class PSO {
public:
    std::vector<Particle> particles;
    std::vector<double> global_best_position;
    std::vector<double> history_global_best_score;
    double global_best_score;
    std::function<double(const std::vector<double>&)> objective_function;  // Objective function stored as a member

    PSO(int swarm_size, int dimensions, std::function<double(const std::vector<double>&)> obj_function,int max_iterations);
    void optimize(int max_iterations, double omega, double phi_p, double phi_g);
    void printResults() const;
};

#endif
