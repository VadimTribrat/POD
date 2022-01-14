#pragma once
#include <vector>
#include <random>
#include <tuple>
#include <functional>
#include "Solver.h"

struct Optimizer
{
	//using func = double(*)(Solver&, std::vector<double>&);
	using func = std::function<double(Solver&, std::vector<double>&)>;
	func f;
	Solver sol;
	std::vector<std::pair<double, double>> constraints;
	std::vector<std::vector<double>> coors;
	unsigned dims;
	double reflect_coeff, compress_coeff, stretch_coeff, global_compress;
	Optimizer(func f, std::vector<std::vector<double>>& init_coors, Solver& sol, const std::vector<std::pair<double, double>>&,
		double rc = 1, double cc = 0.5, double sc = 2, double gc = 0.5);
	std::tuple<size_t, size_t, size_t> points();
	void reflect();
	void compress(const std::vector<double>& x_c, size_t h, size_t l);
	std::vector<double> find_min(int num_iter);
	void fixParams(double&, int);
};
