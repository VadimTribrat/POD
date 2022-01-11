#pragma once
#include <vector>
#include <map>
#include "Point.h"
#include "Lattice.h"
struct Solver
{
	double cohEnergyA;
	std::map<std::pair<Point, Point>, double> distances;
	Lattice lattice;
	Solver(const Lattice& lat, double en);
	std::vector<double> calculateBB(const std::vector<double>& params);
	double energy(const std::vector<double>& params, const std::vector<double>& transform = {});
	double derivative(double energy, const std::vector<double>& params, const std::vector<double>& pos, const std::vector<double>& neg);
};