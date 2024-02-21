#pragma once

#include <set>

double randomDouble ();

bool randomProbability (double probability);

double precision (int low, int high, int numberOfBits);

std::set<int> randomDistinctNumbers (int upperLimit, int amountNumbers);