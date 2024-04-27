#include "differential_evolution.h"
#include "ObjectiveFunctions.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
std::random_device rd;
std::mt19937 random(rd());

DifferentialEvolution::DifferentialEvolution(int population_size, int dimension, int generation_max, double lower_bound, double upper_bound, 
    double F, double Cr)
    : population_size_(population_size), dimension_(dimension), generation_max_(generation_max), lower_bound_(lower_bound)
    , upper_bound_(upper_bound), F_(F), Cr_(Cr) {}

double DifferentialEvolution::runEvolution(const int select, const int recombi) {
    std::vector<double>bestFitnessHistory, worstFitnessHistory, averageFitnessHistory;
    double best=9999999999999999;

    initializePopulation();
    for (size_t g = 0; g < generation_max_; g++){
        
        for (size_t w = 0; w < population_size_; w++){
            
            double fitX=RastriginFunction(population[w].values);
            if(fitX<best){
                best = fitX;
                bestx.values=population[w].values;
                bestx.fitness=fitX;
            }
        }
        
        for(size_t i = 0; i < population_size_;i++){
            std::vector<int> indexs = selection_r();
            int r1=indexs[0];
            int r2=indexs[1];
            int r3=indexs[2];
            
            std::uniform_int_distribution distribution(0, dimension_-1);
            //rand1
            int jrand = distribution(random);
            std::vector<double> child;
            std::uniform_real_distribution<double> distributionR(0.0,1.0);
                
            if(recombi==0){
                for(size_t j = 0; j < dimension_ ;j++){
                //bin
                    double Crrandom = distributionR(random);
                    if (Crrandom<=Cr_ || j==jrand){
                        double x1;
                        if(select==0){
                            x1 = population[r1].values[j];
                        }else if(select==1){
                            x1 = bestx.values[j];
                        }
                        double x2 = population[r2].values[j];
                        double x3 = population[r3].values[j];
                        //se hace la muta con los index
                        auto ui = population[r1].values[j]+F_*(population[r2].values[j]-population[r3].values[j]);
                        child.push_back(ui); 
                    }else{
                        //agregas al hijo el padre
                        child.push_back(population[i].values[j]);
                    }
                }
            }else{
                double Crrandom = distributionR(random);
                int initj=0;
                int j=initj;
                while(Crrandom<=Cr_ && j<dimension_){        
                    double x1;
                    if(select==0){
                        x1 = population[r1].values[j];         
                    }else if(select==1){
                        x1 = bestx.values[j];
                    }
                    double x2 = population[r2].values[j];
                    double x3 = population[r3].values[j];
                    //se hace la muta con los index
                    auto ui = population[r1].values[j]+F_*(population[r2].values[j]-population[r3].values[j]);
                    child.push_back(ui); 
                    distributionR(random);
                    Crrandom = distributionR(random); 
                    j++;
                    initj=j;
                    
                }
                for(size_t j1 = initj; j1 < dimension_ ;j1++){ 
                    child.push_back(population[i].values[j1]);
                }
                
            }

            double fitX=RastriginFunction(population[i].values);
            double fitU=RastriginFunction(child);
            if(fitU<=fitX){
                population[i].values=child;
                population[i].fitness=fitU;
            } 
        }
        if(best==0){
            break;
        }
        
        for (size_t w = 0; w < population_size_; w++){
            double fitX=RastriginFunction(population[w].values);
            if(fitX<best){
                best = fitX;
                bestx.values=population[w].values;
                bestx.fitness=fitX;
            }
        }
        /*
        for(auto x:bestx.values){
            std::cout<<x<<", ";
        }
        std::cout<<std::endl;
        std::cout<<"Fitness: "<<best<<std::endl;
        std::cout<<std::endl;
        */
        
    }
    return bestx.fitness;
    
}

int isRepeat(std::vector<int> index, int num) {
    int boolean=0;
    for(auto i:index){
        if(i==num){
            boolean=1;
            break;
        }
    }
    return boolean;
}

std::vector<int> DifferentialEvolution::selection_r(){
    std::uniform_int_distribution distribution(0, population_size_-1);
    std::vector<int> index;
    for(size_t i=0; i<dimension_;i++){
        int selec=distribution(random);
        if(i!=0){
            int boolean=isRepeat(index,selec);
            if(boolean==1)i--;
            else index.push_back(selec);
        }else{
            index.push_back(selec);
        }     
    }
    return index;
}

void DifferentialEvolution::initializePopulation(){
    for (size_t j = 0; j < dimension_; j++){
        bestx.values.push_back(0.0);
    }
    std::uniform_real_distribution<double> distribution(lower_bound_, upper_bound_);
    for (size_t i = 0; i < population_size_; i++){
        single xi;
        for (size_t j = 0; j < dimension_; j++){
            double dat = distribution(random);
            xi.values.push_back(dat);
        }
        population.push_back(xi);
    }
    
}

