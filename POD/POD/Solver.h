#pragma once
#include <vector>
#include "Lattice.h"

enum class TableParamType {
    a, E_coh, B, C_11, C_12, C_44, E_sol, E_in, E_on, Size
};

struct Solver {
    const double alpha = 1e-6;
    Lattice lattice;
    double energy;
    std::vector<double> target;
    Solver(const Lattice& lattice, const std::vector<double>& target, double energy);

    double Eb(int i, const std::vector<double>& params, const std::vector<double>& transform, bool enableCutoff);
    double Er(int i, const std::vector<double>& params, const std::vector<double>& transform, bool enableCutoff);
    double Energy(const std::vector<double>& params, const std::vector<double>& transform = std::vector<double>(), bool enableCutoff = false);
    std::vector<double> Calculate(const std::vector<double>& params);
    Solver operator=(const Solver& sol);
    Solver() {};
};
