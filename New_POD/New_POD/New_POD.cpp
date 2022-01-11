#include <iostream>
#include <chrono>
#include "Point.h"
#include "Lattice.h"
#include "Solver.h"
#include "Optimizer.h"

double func(Solver& sol, std::vector<double>& params)
{
	double sum_ = 0;
	double param_sum = 0;
	std::vector<double> target{ -2.960, 1.08, 1.32, 0.97, 0.51 };
	std::vector<double> weights{ 0.1, 10, 10, 20, 40 };
	auto res = sol.calculateBB(params);
	for (size_t i = 0; i < res.size(); ++i)
	{
		sum_ += pow(res[i] - target[i], 2) / (target[i] * target[i] * 6);
	}
	//for (size_t i = 0; i < params.size(); ++i)
	//	param_sum += abs(params[i])/1000;
	//for (size_t i = 1; i < params.size(); ++i)
	//{
	//	sum_ += params[i]* params[i] ;
	//}
	//sum_ += (params[0] + 198.9546) * (params[0] + 198.9546);
	return sum_ + param_sum;
}

std::vector<double> GenerateVec(int n, std::vector<std::pair<double, double>>& mid)
{
	std::vector<double> vec(n, 0);

	for (size_t i = 0; i < n; ++i)
	{
		double middle = (mid[i].second + mid[i].first)/2;
		double num = middle + (mid[i].second - mid[i].first)/2 * ((rand() % 10'000) / 5'000.0 - 1.0);
		std::cout << num << " ";
		vec[i] = num;
	}
	std::cout << "\n";
	return vec;

}

int main()
{
	srand(time(0));
	int n = 6;
	//Lattice lattice("Ag.xyz", Point(std::vector<double>{3, 3, 3}), 4.085);
	std::vector<std::vector<double>> init_coors;
	std::vector<std::pair<double, double>> constraints {
		//std::make_pair(0.0, 1.37),
		//std::make_pair(0.0, 0.1),
		//std::make_pair(1.9, 3.3),
		//std::make_pair(7.1, 14),
		//std::make_pair(2.0, 4.0),
		//std::make_pair(0.785, 1.58)
		std::make_pair(0, 1.3),
		std::make_pair(0, 0.3),
		std::make_pair(1.8, 3.3),
		std::make_pair(7.1, 12),
		std::make_pair(-5, 5),
		std::make_pair(0.55, 1.5)
	};
	for (size_t i = 0; i < n+1; ++i)
		init_coors.push_back(GenerateVec(n, constraints));
	Lattice lattice(1,  4.085);
	Solver sol(lattice, 0);
	std::cout << "\nStart:\n";
	auto start = std::chrono::high_resolution_clock::now();

	Optimizer opt(func, init_coors, sol, constraints);
	auto p = opt.find_min(3001);
	auto r = sol.calculateBB(p);

	//auto r = sol.calculateBB({ 0.000856419, 0.107492, 3.44952, 10.9117, 3.9666, 1.46947 });

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::chrono::seconds::period> elapsedtime =
		finish - start;
	std::cout << "elapsed time " << elapsedtime.count() << 's' << std::endl;

	//for (size_t i = 0; i < p.size(); i++)
	//	std::cout << p[i] << " ";
	//std::cout << "\n";

	for (size_t i = 0; i < r.size(); ++i)
		std::cout << r[i] << " ";
}

//0.000856419, 0.107492, 3.44952, 10.9117, 3.9666, 1.46947