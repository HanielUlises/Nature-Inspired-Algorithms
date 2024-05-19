#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "LinearRegression.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::vector<std::vector<double>> readCSV(){
    std::vector<std::vector<double>> DataSet;
    std::ifstream archivo("Salary.csv");
    if (!archivo.is_open()){
        std::cout << "No se pudo abrir el archivo" << std::endl;
    }else{
        std::string linea;
        char delimitador = ',';
        // Leemos la primer línea para descartarla, pues es el encabezado
        getline(archivo, linea);
        std::vector<double> x;
        std::vector<double> y;
        // Leemos todas las líneas
        while (getline(archivo, linea)){

            std::stringstream stream(linea); // Convertir la cadena a un stream
            std::string YearsExperience, Salary;
            // Extraer todos los valores de esa fila
            getline(stream, YearsExperience, delimitador);
            getline(stream, Salary, delimitador);
            x.push_back(std::stod(YearsExperience));
            y.push_back(std::stod(Salary));
        }
        archivo.close();
        DataSet.push_back(x);
        DataSet.push_back(y);
    }
    return DataSet;
}

double evaluation(std::vector<double> ab, std::vector<std::vector<double>> DataSet){
    double score = 0.0f;
    double a = ab[0];
    double b = ab[1];
    std::vector<double> X=DataSet[0];
    std::vector<double> Y_expect=DataSet[1];
    for (size_t i = 0; i < Y_expect.size();++i){
        double y_real = a + b * X[i];
        double dif = y_real - Y_expect[i];
        score+=std::pow(dif,2);
    }
    
    return score;
}

void plottingSLR(std::vector<std::vector<double>> DataSet){
    std::vector<double> X=DataSet[0];
    std::vector<double> Y=DataSet[1];

    std::string title = "Simple Lineal Regression";
    //title.append(function);
    double b = 8555.33918938f;
    double a = 29602.07353482f;

    std::vector<double> Y_fit;
    for (double xi : X) {
        Y_fit.push_back(a + b * xi);
    }    

    plt::scatter(X, Y, 10.0);
    plt::plot(X, Y_fit,{{"color", "red"}});

    plt::xlabel("YearsExperience");
    plt::ylabel("Salary");
    plt::legend();

    plt::title(title);

    plt::show();
    
}

void plotting(const std::vector<std::vector<double>>& history) {
    if (history.empty()) {
        std::cerr << "Error: History is empty. No data to plot.\n";
        return;
    }

    size_t history_size = history[0].size();
    std::vector<double> bestscore = history[0];
    std::vector<double> worstscore = history[0];
    for (const auto &score : history) {
        if (score.size() != history_size) {
            std::cerr << "Error: Inconsistent sizes in history data.\n";
            return;
        }
        if (score.back() <= bestscore.back()) bestscore = score;
        if (score.back() >= worstscore.back()) worstscore = score;
    }

    std::string title = "Convergence Graph";

    std::vector<int> generations(history_size);
    std::iota(generations.begin(), generations.end(), 0);

    plt::plot(generations, bestscore, {{"label", "best"}});
    plt::plot(generations, worstscore, {{"label", "worst"}});

    plt::xlabel("Generation");
    plt::ylabel("Best global score");
    plt::legend();

    plt::title(title);
    plt::show();
}


void plottingSLR_withSolution(std::vector<std::vector<double>> DataSet, std::vector<double> particle){
    std::vector<double> X=DataSet[0];
    std::vector<double> Y=DataSet[1];

    std::string title = "Simple Lineal Regression";
    //title.append(function);
    double b = 8555.33918938f;
    double a = 29602.07353482f;

    std::vector<double> Y_fit;
    for (double xi : X) {
        Y_fit.push_back(a + b * xi);
    }    

    double a_particle = particle[0];
    double b_particle = particle[1];

    std::vector<double> Y_minimize;
    for (double xi : X) {
        Y_minimize.push_back(a_particle + b_particle * xi);
    }  

    plt::scatter(X, Y, 10.0);
    plt::plot(X, Y_fit,{{"color", "red"}});
    plt::plot(X, Y_minimize,{{"color", "blue"}});

    plt::xlabel("YearsExperience");
    plt::ylabel("Salary");
    plt::legend();

    plt::title(title);

    plt::show();
    
}