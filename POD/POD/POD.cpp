#include <iostream>
#include <omp.h>
#include <cmath>
#include <functional>
#include <random>
#include <chrono>
#include "Lattice.h"
#include "Solver.h"
#include "Optimizer.h"


double func(Solver& sol, std::vector<double>& params)
{
	double sum_ = 0;
	std::vector<double> target{-2.960, 1.08, 1.32, 0.97, 0.51};
	std::vector<double> weights{ 0.1, 10, 10, 20, 40};
	auto res = sol.Calculate(params);
	for (size_t i = 0; i < res.size(); ++i)
	{
		if (isnan(res[i]))
		{
			std::cout << "\n\n\nAlarm: " << i << "\n\n\n";
		}
		sum_ += pow(res[i] - target[i], 2)/(target[i]*target[i]);
	} 
	//for (size_t i = 1; i < params.size(); ++i)
	//{
	//	sum_ += params[i]* params[i] ;
	//}
	//sum_ += (params[0] + 198.9546) * (params[0] + 198.9546);
	return sum_;
}
std::vector<double> GenerateVec(int n)
{
	std::vector<double> vec(n, 0);
	std::cout << "Vector init: \n";
	
	for (size_t i = 0; i < n; ++i)
	{
		double num = (rand() % 2'000) /2000.0;
		std::cout << num << " ";
		vec[i] = num;
	}
	std::cout << "\n";
		return vec;

}

int main()
{
	srand(time(0));
	Lattice lattice("Ag.xyz", Point(std::vector<double>{3, 3, 3}), 4.085);
	std::cout << "Periodic vector " << lattice.period << "\n";
	Solver sol(lattice, std::vector<double>{0}, 0.5);
	std::vector<std::vector<double>> init_;

	//auto p = std::vector<double>{ 10.8772, -0.91834, -0.0086409, 2.68355, 3.89788, 18.9771 };
	for (size_t i = 0; i < 7; ++i)
	{
		init_.push_back(GenerateVec(6));
	}

	std::cout << "\nStart:\n";
	auto start = std::chrono::high_resolution_clock::now();
	Optimizer opt(func, init_, sol);
	auto p = opt.find_min(101);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::chrono::seconds::period> elapsedTime =
		finish - start;
	std::cout << "Elapsed time " << elapsedTime.count() << 's' << std::endl;


	//auto r = sol.Calculate(p);
	//std::cout << "\n\n" << func(sol, p) << "\n\n";
	//for (size_t i = 0; i < p.size(); ++i)
	//	std::cout << p[i] << " ";
	//std::cout << "\n";

	//for (size_t i = 0; i < r.size(); ++i)
	//	std::cout << r[i] << " ";
	//std::cout << "\n";

	return 0;
}
