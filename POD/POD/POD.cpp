#include <iostream>
#include <omp.h>
#include <cmath>
#include <functional>
#include "Lattice.h"
#include "Solver.h"
#include "Optimizer.h"

double func(Solver& sol, std::vector<double>& params)
{
	return std::abs(sol.Calculate(params)[0] + 2.960);
}

int main()
{
	Lattice lattice("Ag.xyz", Point(std::vector<double>{4.085, 4.085, 4.085}), 1);
	Solver sol(lattice, std::vector<double>{0}, 0.5);
	std::vector<std::vector<double>> init_;
	std::cout << sol.Energy(std::vector<double>{0.2, 0.1523, 0.405, 0.678, 0.32, 0.1}) << "\n";
	auto p = std::vector<double>{ 0.277627, 0.00595991, 0.47706, 0.683957, 0.33074, 0.259035 };
	/*std::cout << func(solver, p);*/
	auto r =  sol.Calculate(std::vector<double>{0.18953, 1.999, 4.12345, 0.6436, 0.32769, 0.1296});
	for (size_t i = 0; i < r.size(); ++i)
		std::cout << r[i] << " ";
	std::cout << "\n";
	init_.push_back(std::vector<double>{0.2, 0.1523, 0.405, 0.675, 0.321435, 0.175647});
	init_.push_back(std::vector<double>{0.2549352, 0.18694, 0.405, 0.6145, 0.332345, 0.1999});
	init_.push_back(std::vector<double>{0.21234, 0.1235, 0.465, 0.678, 0.3198435, 0.13443});
	init_.push_back(std::vector<double>{0.18953, 0.1999, 0.412345, 0.6436, 0.32769, 0.1296});
	init_.push_back(std::vector<double>{0.123753, 0.11345, 0.39995, 0.63457, 0.329999, 0.1543869});
	init_.push_back(std::vector<double>{0.259, 0.187462, 0.4356, 0.6666, 0.37854, 0.15386972});
	init_.push_back(std::vector<double>{0.24509, 0.1555, 0.38678, 0.67777, 0.36549387, 0.12478629});

	for (auto& val : init_)
	{
		auto r = sol.Calculate(val);
		for (size_t i = 0; i < r.size(); ++i)
			std::cout << r[i] << " ";
		std::cout << "\n";
	}

	Optimizer opt(func, init_, sol);
	opt.find_min(10001);
	std::cout << sol.Calculate(p)[0];

	return 0;
}
