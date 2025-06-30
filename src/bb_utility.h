#pragma once
#include <vector>

int combination(int n, int k);
double parallel_reliability(int n, double r);
double parallel_reliability(const std::vector<double>& reliabilities);
double k_out_of_n_reliability(int n, int k, double r);