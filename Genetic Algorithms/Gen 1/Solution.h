#pragma once

#include <vector>
#include <string>

class Solution{
    public:
    Solution(int numberOfBits, int low, int high);
    std::string toString();
    double bitsToDouble ();

    private:
    int mNumberOfBits, mLow, mHigh;
    std::vector<int> bits; // vector of chromosomes

};