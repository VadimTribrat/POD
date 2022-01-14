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
double energy(const std::vector<double>& params, double dist);

int main()
{
	srand(time(0));
	int n = 6;
	//Lattice lattice("Ag.xyz", Point(3, 3, 3), 4.085);
	std::vector<std::vector<double>> init_coors;
	std::vector<std::pair<double, double>> constraints{
		std::make_pair(-0.2, 0.1),
		std::make_pair(0, 0.2),
		std::make_pair(1.9, 3.0),
		std::make_pair(10, 12),
		std::make_pair(2, 5),
		std::make_pair(0.9, 1.3)
		//std::make_pair(-0.3, 0.1),
		//std::make_pair(0, 0.2),
		//std::make_pair(1.9, 3.3),
		//std::make_pair(9, 13),
		//std::make_pair(2, 5),
		//std::make_pair(0.9, 2)
	};
	std::vector<std::pair<double, double>> constraints1
	{
		std::make_pair(-0.3, 0.1),
		std::make_pair(0, 0.2),
		std::make_pair(1.9, 3.3),
		std::make_pair(9, 13),
		std::make_pair(2, 5),
		std::make_pair(0.9, 2)
	};

	for (size_t i = 0; i < n+1; ++i)
		init_coors.push_back(GenerateVec(n, constraints));
	Lattice lattice(Point(3, 3, 3),  4.085);
	Solver sol(lattice, -4.32);


	std::cout << "\nStart:\n";
	auto start = std::chrono::high_resolution_clock::now();
	auto f = std::bind(func, _1, _2, std::vector<double>{ -2.960, 1.08, 1.32, 0.97, 0.51 }, &Solver::calculateBB);
	Optimizer opt(f, init_coors, sol, constraints);
	auto p = opt.find_min(251);
	auto r = sol.calculateBB(p);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
		finish - start;
	std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	printVec(p);
	std::cout << sol.lattice.mul << " ";
	printVec(r);

	std::cout << "\nStart:\n";
	start = std::chrono::high_resolution_clock::now();
	sol.E_cohB = r[0];
	sol.lattice.setMul(sol.lattice.mul);
	sol.fittedBB = p;
	std::cout << " " << sol.lattice.mul << " ";
	Optimizer opt1(std::bind(func, _1, _2, std::vector<double>{0.63}, & Solver::calculateAB), init_coors, sol, constraints);
	auto p1 = opt1.find_min(501);
	auto r1 = sol.calculateAB(p1);
	finish = std::chrono::high_resolution_clock::now();
	elapsedtime =
		finish - start;
	std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	printVec(p1);
	printVec(r1);
	sol.fittedBB = std::vector<double>{};

	std::cout << "\nStart:\n";
	start = std::chrono::high_resolution_clock::now();
	sol.fittedBB = p;
	sol.fittedAB = p1;
	sol.E_cohB = r[0];
	sol.lattice.setMul(sol.lattice.mul);
	Optimizer opt2(std::bind(func, _1, _2, std::vector<double>{ -0.08, -0.69 }, & Solver::calculateAA), init_coors, sol, constraints1);
	auto p2 = opt2.find_min(501);
	auto r2 = sol.calculateAA(p2);
	finish = std::chrono::high_resolution_clock::now();
	elapsedtime =
		finish - start;
	std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	sol.fittedBB = std::vector<double>{};
	sol.fittedAB = std::vector<double>{};
	printVec(p2);
	printVec(r2);

	//std::cout << "\nStart:\n";
	//auto start = std::chrono::high_resolution_clock::now();
	//sol.E_cohB = -2.96104;
	// sol.lattice.setMul(4.02547);
	//sol.fittedBB = std::vector<double>{ -0.0159319, 0.111505, 2.82208, 10.2265, 3.0267, 1.20272 };
	//std::cout << " " << sol.lattice.mul << " ";
	//Optimizer opt(std::bind(func, _1, _2, std::vector<double>{0.63}, &Solver::calculateAB), init_coors, sol, constraints);
	//auto p = opt.find_min(101);
	//auto r = sol.calculateAB(p);
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
	//	finish - start;
	//std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	//printVec(p);
	//printVec(r);
	//sol.fittedBB = std::vector<double>{};

	//std::cout << "\nStart:\n";
	//auto start = std::chrono::high_resolution_clock::now();
	//sol.E_cohB = -2.96104;
	//sol.lattice.setMul(4.02547);
	//sol.fittedBB = std::vector<double>{ -0.0159319, 0.111505, 2.82208, 10.2265, 3.0267, 1.20272 };
	//sol.fittedAB = std::vector<double>{ -0.057879, 0.102937, 2.69373, 10.3437, 2.36513, 1.23069, 0.629993 };
	//Optimizer opt = Optimizer(std::bind(func, _1, _2, std::vector<double>{ -0.08, -0.69 }, &Solver::calculateAA), init_coors, sol, constraints);
	//auto p = opt.find_min(301);
	//auto r = sol.calculateAA(p);
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
	//	finish - start;
	//std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;
	//sol.fittedBB = std::vector<double>{};
	//sol.fittedAB = std::vector<double>{};
	//printVec(p);
	//printVec(r);

	//sol.fittedBB = std::vector<double>{ -0.0159319, 0.111505, 2.82208, 10.2265, 3.0267, 1.20272 };
	//sol.lattice.setMul(4.02547);
	//printVec(sol.calculateBB({}));

	//sol.E_cohB = -2.96104;
	//sol.lattice.setMul(4.02547);
	//sol.fittedBB = std::vector<double>{ -0.0159319, 0.111505, 2.82208, 10.2265, 3.0267, 1.20272 };
	//printVec(sol.calculateAB({ -0.057879, 0.102937, 2.69373, 10.3437, 2.36513, 1.23069, 0.629993 }));
	//sol.fittedBB = std::vector<double>{};

	//sol.E_cohB = -2.96104;
	//sol.lattice.setMul(4.02547);
	//sol.fittedBB = std::vector<double>{ -0.0159319, 0.111505, 2.82208, 10.2265, 3.0267, 1.20272 };
	//sol.fittedAB = std::vector<double>{ -0.057879, 0.102937, 2.69373, 10.3437, 2.36513, 1.23069, 0.629993 };
	//printVec(sol.calculateAA({ -0.0580274, 0.060168, 3.29989, 10.0452, 2.93416, 1.25747 }));
	//sol.fittedBB = std::vector<double>{};
	//sol.fittedAB = std::vector<double>{};
}


//Params BB: -0.0159319, 0.111505, 2.82208, 10.2265, 3.0267, 1.20272
//Params AB: -0.057879, 0.102937, 2.69373, 10.3437, 2.36513, 1.23069
//Params AA: -0.0580274, 0.060168, 3.29989, 10.0452, 2.93416, 1.25747

// BB_target: 0.599535, 0.289139, 2.71249, 11.8106, 0.532275, 0.881892
// AB_target: 0.63
// AA_target: -0.08, -0.69

double func(Solver& sol, std::vector<double>& params, const std::vector<double> target, std::vector<double> (Solver::*func)(const std::vector<double>& ))
{
	double sum_ = 0;
	auto res = (sol.*func)(params);
	for (size_t i = 0; i < res.size(); ++i)
	{
		sum_ += pow(res[i]/target[i] - 1, 2) / (target.size());
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

double energy(const std::vector<double>& params, double distance)
{
	return ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1))) -
		sqrt(params[5] * params[5] * exp(-2 * params[4] * (distance / params[2] - 1)));
}