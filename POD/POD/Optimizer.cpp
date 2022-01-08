#include "Optimizer.h"
#include <iostream>

Optimizer::Optimizer(func f, std::vector<std::vector<double>>& init_coors, Solver& sol,
	double rc, double cc, double sc, double gc): sol(sol)
{
	this->f = f;
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
	}
	for (size_t i = 0; i < dims; ++i)
	{
		x_r[i] = (1 + reflect_coeff) * x_c[i] - reflect_coeff * coors[h][i];
	}
	//std::cout << "x_r: \n";
	//for (size_t i = 0; i < dims; ++i)
	//	std::cout << x_r[i] << " ";
	//std::cout << "\n";
	//std::cout << "f(sol, x_r) = " << f(sol, x_r) << "\n";
	if (f(sol, x_r) < f(sol, coors[l]))
	{
		std::vector<double> x_e(dims, 0);
		for (size_t i = 0; i < size - 1; ++i)
		{
			x_e[i] = (1 - stretch_coeff) * x_c[i] + stretch_coeff * x_r[i];
		}
		if (f(sol, x_e) < f(sol, x_r))
			coors[h] = x_e;
		else
			coors[h] = x_r;
	}
	else if ((f(sol, x_r) >= f(sol, coors[l])) and (f(sol, x_r) < f(sol, coors[g])))
	{
		coors[h] = x_r;
	}
	else if ((f(sol, x_r) >= f(sol, coors[g])) and (f(sol, x_r) < f(sol, coors[h])))
	{
		coors[h] = x_r;
		compress(x_c, h, l);
	}
	else if (f(sol, coors[h]) <= f(sol, x_r))
	{
		compress(x_c, h, l);
	}
	else
		std::cout << "Compress\n";
}

void Optimizer::compress(const std::vector<double> x_c, size_t h, size_t l)
{
	std::vector<double> x_s(dims, 0);
	auto x_h = coors[h];
	for (size_t i = 0; i < dims; ++i)
	{
		x_s[i] = (1 - compress_coeff) * x_c[i] + compress_coeff * x_h[i];
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
			}
		}
	}
}

std::vector<double> Optimizer::find_min(int num_iter)
{
	for (size_t i = 0; i < num_iter; ++i)
	{
		reflect();
		auto res = coors[std::get<0>(points())];
		if (i % 10 == 0)
		{
			for (size_t i = 0; i < res.size(); ++i)
				std::cout << res[i] << " ";
			std::cout << "\n";
		}
	}
	return coors[std::get<0>(points())];
}