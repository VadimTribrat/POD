#pragma once
#include <vector>
#include <map>
#include "Point.h"
#include "Lattice.h"
struct Solver
{
	std::vector<double> fittedBB;
	std::vector<double> fittedAB;
	std::vector<double> fittedAA;
	Lattice lattice;
	double E_cohA, E_cohB;
	Solver(const Lattice& lat, double e);
	std::vector<double> calculateBB(const std::vector<double>& params);
	std::vector<double> calculateAB(const std::vector<double>& params);
	std::vector<double> calculateAA(const std::vector<double>& params);
	double energy(const std::vector<double>& params, const std::vector<double>& transform = {});
	double derivative(double energy, const std::vector<double>& params, const std::vector<double>& pos, const std::vector<double>& neg);
};