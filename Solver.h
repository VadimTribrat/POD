#pragma once
#include <vector>
#include "Lattice.h"

struct Solver {
    const double alpha = 0.1;
    Lattice lattice;
    double energy;
    std::vector<double> target;
    Solver(const Lattice& lattice, const std::vector<double>& target, double energy);

    double Eb(int i, const std::vector<double>& params, const std::vector<double>& transform, bool enableCutoff);
    double Er(int i, const std::vector<double>& params, const std::vector<double>& transform, bool enableCutoff);
    double Energy(const std::vector<double>& params,
        const std::vector<double>& transform = std::vector<double>(), bool enableCutoff = true);
    std::vector<double> Calculate(const std::vector<double>& params);
    Solver operator=(const Solver& sol);
    Solver() {};
    double Derivative2(double E, const std::vector<double>& params, const std::vector<double>& positiveTransform,
        const std::vector<double>& negativeTransform);
};
