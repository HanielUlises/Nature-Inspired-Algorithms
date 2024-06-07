#!/bin/bash

g++ -o aco_algorithm main.cpp ACO.cpp -std=c++11

./aco_algorithm

if [ -s ACO_Results.csv ]; then
    python3 ANOVA.py
else
    echo "El archivo CSV no se generó correctamente."
fi