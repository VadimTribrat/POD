#include "Solver.h"
#include <cmath>
#include <omp.h>

Solver::Solver(const Lattice& lattice, const std::vector<double>& target, double energy) :
	lattice(lattice), target(target), energy(energy)
{
}

double Solver::Er(int i, const std::vector<double>& params, const std::vector<double>& transform, bool enableCutoff)
{
	double sum_ = 0;
	for (size_t j = 0; j < lattice.atoms.size(); ++j)
	{
		if (i == j)
			continue;
		double distance = lattice.dist(i, j, transform, enableCutoff);
		sum_ += ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1)));
	}
	return sum_;
}
double Solver::Eb(int i, const std::vector<double>& params, const std::vector<double>& transform, bool enableCutoff) 
{
	double sum_ = 0;
	for (size_t j = 0; j < lattice.atoms.size(); ++j)
	{
		if (i == j)
			continue;
		double distance = lattice.dist(i, j, transform, enableCutoff);
		sum_ += -sqrt(params[5]*params[5]*exp(-2*params[4]*(distance/params[2] - 1)));
	}
	return sum_;
}

double Solver::Energy(const std::vector<double>& params,
	const std::vector<double>& transform, bool enableCutoff)
{
	double sum_ = 0;
	for (size_t i = 0; i < lattice.atoms.size(); ++i)
	{
		//std::cout << "i = " << i <<" " << Eb(i, params, std::vector<double>(), enableCutoff)
		//	+ Er(i, params, std::vector<double>(), enableCutoff) << "\n";
		sum_ += Eb(i, params, std::vector<double>(), enableCutoff) + Er(i, params, std::vector<double>(), enableCutoff);
	}
	return sum_;

}

std::vector<double> Solver::Calculate(const std::vector<double>& params)
{
	std::vector<double> res;
	res.push_back(Energy(params)/lattice.atoms.size());
	return res;
}

Solver Solver::operator=(const Solver& sol)
{
	Solver temp;
	temp.lattice = sol.lattice;
	temp.energy = sol.energy;
	temp. target = sol.target;
	return temp;
}