#include "Optimizer.h"
#include <iostream>

Optimizer::Optimizer(func f, std::vector<std::vector<double>>& init_coors,
	Solver& sol, const std::vector<std::pair<double, double>>& cons, 
	double rc, double cc, double sc, double gc) : sol(sol), constraints(cons), f(f)
{
	reflect_coeff = rc;
	compress_coeff = cc;
	stretch_coeff = sc;
	global_compress = gc;
	coors = init_coors;
	dims = init_coors[0].size();
}

std::tuple<size_t, size_t, size_t> Optimizer::points()
{
	std::vector<std::pair<int, double>> values(dims + 1, { 0, 0 });
	size_t h, g, l;
	for (size_t i = 0; i < coors.size(); ++i)
		values[i] = { i, f(sol, coors[i]) };
	std::stable_sort(values.begin(), values.end(),
		[](const std::pair<int, double >& a, const std::pair<int, double >& b)
		{return a.second > b.second; });
	h = values[0].first;
	g = values[1].first;
	l = values[values.size() - 1].first;
	return std::make_tuple(h, g, l);
}

void Optimizer::reflect()
{
	size_t size = dims + 1, h, g, l;
	auto res = points();
	std::vector<double> x_c(dims, 0);
	std::vector<double> x_r(dims, 0);
	h = std::get<0>(res), g = std::get<1>(res), l = std::get<2>(res);
	for (size_t i = 0; i < size; i++)
		for (size_t j = 0; j < size - 1; j++)
			x_c[j] += coors[i][j];
	for (size_t i = 0; i < size - 1; ++i)
	{
		x_c[i] -= coors[h][i];
		x_c[i] /= size - 1;
		fixParams(x_c[i], i);
	}
	for (size_t i = 0; i < dims; ++i)
	{
		x_r[i] = (1 + reflect_coeff) * x_c[i] - reflect_coeff * coors[h][i];
		fixParams(x_r[i], i);
	}
	auto result = f(sol, x_r);
	if (result < f(sol, coors[l]))
	{
		std::vector<double> x_e(dims, 0);
		for (size_t i = 0; i < size - 1; ++i)
		{
			x_e[i] = (1 - stretch_coeff) * x_c[i] + stretch_coeff * x_r[i];
			fixParams(x_e[i], i);
		}
		if (f(sol, x_e) < f(sol, x_r))
			coors[h] = x_e;
		else
			coors[h] = x_r;
	}
	else if ((result >= f(sol, coors[l])) and (result < f(sol, coors[g])))
		coors[h] = x_r;
	else if ((result >= f(sol, coors[g])) and (result < f(sol, coors[h])))
	{
		coors[h] = x_r;
		compress(x_c, h, l);
	}
	else if (f(sol, coors[h]) <= result)
		compress(x_c, h, l);
	else
	{
		compress(x_c, h, l);
	}
}

void Optimizer::compress(const std::vector<double> x_c, size_t h, size_t l)
{
	//std::cout << "\n\nCompress\n\n";
	std::vector<double> x_s(dims, 0);
	auto x_h = coors[h];
	for (size_t i = 0; i < dims; ++i)
	{
		x_s[i] = (1 - compress_coeff) * x_c[i] + compress_coeff * x_h[i];
		fixParams(x_s[i], i);
	}
	if (f(sol, x_s) < f(sol, x_h))
		coors[h] = x_s;
	else
	{
		for (size_t i = 0; i < dims + 1; ++i)
		{
			for (size_t j = 0; j < dims; ++j)
			{
				coors[i][j] = coors[l][j] + global_compress * (coors[i][j] - coors[l][j]);
				fixParams(coors[i][j], j);
			}
		}
	}
}

std::vector<double> Optimizer::find_min(int num_iter)
{
	for (size_t i = 0; i < num_iter; ++i)
	{
		reflect();
		if (i % 100 == 0)
		{
			std::cout << "Iter:\n";
			auto poin = coors[std::get<0>(points())];
			auto res = sol.calculateBB(poin);
			for (size_t i = 0; i < poin.size(); ++i)
				std::cout << poin[i] << ", ";
			std::cout << "\n";
			for (size_t i = 0; i < res.size(); ++i)
				std::cout << res[i] << ", ";
			std::cout << "\n" << f(sol, poin) << "\n\n";
		}
		for (auto& val : coors)
			for (size_t i = 0; i < val.size(); ++i)
				fixParams(val[i], i);
	}
	return coors[std::get<0>(points())];
}

void Optimizer::fixParams(double& x, int i)
{

	if (x <= constraints[i].first || x >= constraints[i].second)
	{
		if (x <= constraints[i].first)
			x = constraints[i].first + (constraints[i].second - constraints[i].first) / 5 * ((rand() % 10'000) / 10'000.0);
		if (x >= constraints[i].second)
			x = constraints[i].second - (constraints[i].second - constraints[i].first) / 5 * ((rand() % 10'000) / 10'000.0);
	}
}
