#include <cmath>
#include "bb_utility.h"

int combination(const int n, const int k) {
    int result = 1;
    if (k > n) return 0;

    for (int i = 1; i <= k; ++i) {
        result = result * (n - i + 1) / i;
    }

    return result;
}

double parallel_reliability(const int n, const double r) {
    return 1.0 - std::pow(1 - r, n);
}

double parallel_reliability(const std::vector<double> &reliabilities) {
    double failure_probability = 1.0;

    for (double r : reliabilities) {
        failure_probability *= (1 - r);
    }

    return 1.0 - failure_probability;
}

double k_out_of_n_reliability(const int n, const int k, const double r) {
    double reliability = 0.0;
    for (int i = k; i <= n; ++i) {
        auto coeff = combination(n, i);
        reliability += static_cast<double>(coeff) * std::pow(r, i) * std::pow(1.0 - r, n - i);
    }
    return reliability;
}
