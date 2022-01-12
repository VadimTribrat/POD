#include <iostream>
#include <chrono>
#include <functional>
#include "Point.h"
#include "Lattice.h"
#include "Solver.h"
#include "Optimizer.h"

using namespace std::placeholders;

double func(Solver& sol, std::vector<double>& params, const std::vector<double> target, std::vector<double>(Solver::* func)(const std::vector<double>&));
std::vector<double> GenerateVec(int n, std::vector<std::pair<double, double>>& mid);
void printVec(const std::vector<double>&);

int main()
{
	srand(time(0));
	int n = 6;
	//Lattice lattice("Ag.xyz", Point(3, 3, 3), 4.085);
	std::vector<std::vector<double>> init_coors;
	std::vector<std::pair<double, double>> constraints {
		std::make_pair(-2, 1.3),
		std::make_pair(0, 0.3),
		std::make_pair(1.8, 3.3),
		std::make_pair(7.1, 12),
		std::make_pair(-5, 5),
		std::make_pair(0.55, 1.5)
	};
	for (size_t i = 0; i < n+1; ++i)
		init_coors.push_back(GenerateVec(n, constraints));
	Lattice lattice(Point(3, 3, 3),  4.085);
	Solver sol(lattice, -4.32);


	//std::cout << "\nStart:\n";
	//auto start = std::chrono::high_resolution_clock::now();
	//auto f = std::bind(func, _1, _2, std::vector<double>{ -2.960, 1.08, 1.32, 0.97, 0.51 }, &Solver::calculateBB);
	//Optimizer opt(f, init_coors, sol, constraints);
	//auto p = opt.find_min(201);
	//auto r = sol.calculateBB(p);
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
	//	finish - start;
	//std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	//printVec(p);
	//std::cout << sol.lattice.multiplier << " ";
	//printVec(r);

	//printVec(sol.calculateBB({ -0.184678, 0.132579, 2.84169, 8.13273, 3.22019, 1.24989 }));



	std::cout << "\nStart:\n";
	auto start = std::chrono::high_resolution_clock::now();
	sol.E_cohB = -2.95;
	sol.fittedBB = std::vector<double>{ -0.184678, 0.132579, 2.84169, 8.13273, 3.22019, 1.24989 };
	Optimizer opt(std::bind(func, _1, _2, std::vector<double>{0.63}, &Solver::calculateAB), init_coors, sol, constraints);
	auto p = opt.find_min(201);
	auto r = sol.calculateAB(p);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
		finish - start;
	std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	printVec(p);
	printVec(r);
	sol.fittedBB = std::vector<double>{};



	//std::cout << "\nStart:\n";
	//auto start = std::chrono::high_resolution_clock::now();
	//sol.fittedBB = std::vector<double>{ -0.184678, 0.132579, 2.84169, 8.13273, 3.22019, 1.24989 };
	//sol.fittedAB = std::vector<double>{ 0.0364029, 0.162998, 3.07157, 7.79454, 3.06882, 1.43693 };
	//Optimizer opt = Optimizer(std::bind(func, _1, _2, std::vector<double>{ -0.08, -0.69 }, &Solver::calculateAA), init_coors, sol, constraints);
	//auto p = opt.find_min(301);
	//auto r = sol.calculateAA(p);
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
	//	finish - start;
	//std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	////auto p = std::vector<double>{ -0.310013, 0.0875706, 2.62258, 9.27129, 3.78214, 1.13736 };
	//sol.fittedBB = std::vector<double>{};
	//sol.fittedAB = std::vector<double>{};
	//printVec(p);
	//printVec(r);
}


//Params BB: -0.184678, 0.132579, 2.84169, 8.13273, 3.22019, 1.24989
//Params AB: 0.0364029, 0.162998, 3.07157, 7.79454, 3.06882, 1.43693
//Params AA: -1.9912 0.256259 2.54045 8.86415 0.593601 1.49264

// BB_target: 0.599535, 0.289139, 2.71249, 11.8106, 0.532275, 0.881892
// AB_target: 0.63
// AA_target: -0.08, -0.69

double func(Solver& sol, std::vector<double>& params, const std::vector<double> target, std::vector<double> (Solver::*func)(const std::vector<double>& ))
{
	double sum_ = 0;
	auto res = (sol.*func)(params);
	for (size_t i = 0; i < res.size(); ++i)
	{
		sum_ += pow(res[i] - target[i], 2) / (target[i] * target[i] * target.size());
	}
	return sum_;
}

std::vector<double> GenerateVec(int n, std::vector<std::pair<double, double>>& mid)
{
	std::vector<double> vec(n, 0);
	for (size_t i = 0; i < n; ++i)
	{
		double middle = (mid[i].second + mid[i].first) / 2;
		double num = middle + (mid[i].second - mid[i].first) / 2 * ((rand() % 10'000) / 5'000.0 - 1.0);
		vec[i] = num;
	}
	return vec;
}

void printVec(const std::vector<double>& vec)
{
	for (size_t i = 0; i < vec.size(); ++i)
		std::cout << vec[i] << " ";
	std::cout << "\n";
}