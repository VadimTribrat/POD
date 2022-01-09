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
		if (distance != 0)
		{
			sum_ += ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1)));
		}
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
		if (distance != 0)
			sum_ += params[5] * params[5] * exp(-2 * params[4] * (distance / params[2] - 1));
	}
	return -sqrt(sum_);
}

double Solver::Energy(const std::vector<double>& params,
	const std::vector<double>& transform, bool enableCutoff)
{
	double sum_ = 0;
#pragma omp parallel num_threads(8)
	{
#pragma omp for
		for (int i = 0; i < lattice.atoms.size(); ++i)
		{
#pragma omp atomic
			sum_ += Eb(i, params, transform, enableCutoff) + Er(i, params, transform, enableCutoff);
		}
	}
	if (isnan(sum_))
	{
		std::cout << "Energy  nan\n\n";
	}
	return sum_;
}

std::vector<double> Solver::Calculate(const std::vector<double>& params)
{
	std::vector<double> res;
	double E_0 = Energy(params, {}, true);
	double E_coh = E_0 / lattice.atoms.size();
	//double E_coh = -2.960;
	if (isnan(E_coh))
	{
		std::cout << "E_coh  nan\n\n";
		E_coh = std::numeric_limits<double>::max();
	}

	double V_0 = std::pow(lattice.multiplier, 3) * 3 * 3 * 3;

	double C_11_derivative = Derivative2(E_0, params, { 1 + alpha, 1 + alpha, 1 }, { 1 - alpha, 1 - alpha, 1 });
	double C_12_derivative = Derivative2(E_0, params, { 1 + alpha, 1 - alpha, 1 }, { 1 - alpha, 1 + alpha, 1 });

	double C_11 = 1.602*((C_11_derivative + C_12_derivative) / (4 * V_0));
	double C_12 = 1.602*((C_11_derivative - C_12_derivative) / (4 * V_0));

	double B_derivative = Derivative2(E_0, params, { 1 + alpha, 1 + alpha, 1 + alpha }, { 1 - alpha, 1 - alpha, 1 - alpha });
	double B = 1.602*(B_derivative / (9 * V_0));

	double C_44_derivative = Derivative2(E_0, params, { 1, alpha, alpha, 1, 1 / (1 - alpha * alpha) },
		{ 1, -alpha, -alpha, 1, 1 / (1 - alpha * alpha) });
	double C_44 = 1.602*(C_44_derivative / (4 * V_0));
	if (isnan(B))
	{
		std::cout << "B  nan\n\n";
		B = std::numeric_limits<double>::max();
	}
	if (isnan(C_11))
	{
		std::cout << "C_11  nan\n\n";
		C_11 = std::numeric_limits<double>::max();
	}
	if (isnan(C_12))
	{
		std::cout << "C_12  nan\n\n";
		C_12 = std::numeric_limits<double>::max();
	}
	if (isnan(C_44))
	{
		std::cout << "C_44  nan\n\n";
		C_44 = std::numeric_limits<double>::max();
	}
	res.push_back(E_coh);
	res.push_back(B);
	res.push_back(C_11);
	res.push_back(C_12);
	res.push_back(C_44);
	return res;
}

Solver Solver::operator=(const Solver& sol)
{
	Solver temp;
	temp.lattice = sol.lattice;
	temp.energy = sol.energy;
	temp.target = sol.target;
	return temp;
}

double Solver::Derivative2(double E, const std::vector<double>& params, const std::vector<double>& positiveTransform,
	const std::vector<double>& negativeTransform) {
	double positiveEnergy = Energy(params, positiveTransform);
	double negativeEnergy = Energy(params, negativeTransform);

	return (positiveEnergy - 2 * E + negativeEnergy) / (alpha * alpha);
}