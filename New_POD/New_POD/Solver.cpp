#include <omp.h>
#include "Solver.h"

Solver::Solver(const Lattice& lat, double en):
	lattice(lat), cohEnergyA(en)
{
	//for (size_t i = 0; i < lattice.atoms.size(); ++i)
	//	for (size_t j = i + 1; j < lattice.atoms.size(); ++j)
	//		distances[std::make_pair(lattice[i].second, lattice[j].second)] = lattice.distance(i, j);
}


std::vector<double> Solver::calculateBB(const std::vector<double>& params)
{
	std::vector<double> res;
	double E_0 = energy(params);
	double E_coh = E_0 / lattice.atoms.size();

	double V_0 = std::pow(lattice.multiplier, 3) * lattice.period[0] * lattice.period[1] * lattice.period[2];

	double C_11_derivative = derivative(E_0, params, { 1 + 1e-2, 1 + 1e-2, 1 }, { 1 - 1e-2, 1 - 1e-2, 1 });
	double C_12_derivative = derivative(E_0, params, { 1 + 1e-2, 1 - 1e-2, 1 }, { 1 - 1e-2, 1 + 1e-2, 1 });

	double C_11 = ((C_11_derivative + C_12_derivative) / (2*V_0))*0.8019;
	double C_12 = ((C_11_derivative - C_12_derivative) / (2 * V_0))*0.8019;

	double B_derivative = derivative(E_0, params, { 1 + 1e-2, 1 + 1e-2, 1 + 1e-2 }, { 1 - 1e-2, 1 - 1e-2, 1 - 1e-2 });
	double B = 2*(B_derivative / (9 * V_0))*0.8019;

	double C_44_derivative = derivative(E_0, params, { 1, 1e-2, 1e-2, 1, 1 / (1 - 1e-2 * 1e-2) },
		{ 1, -1e-2, -1e-2, 1, 1 / (1 - 1e-2 * 1e-2) });
	double C_44 = (C_44_derivative / (2 * V_0))*0.8019;
	res.push_back(E_coh);
	res.push_back(B);
	res.push_back(C_11);
	res.push_back(C_12);
	res.push_back(C_44);
	return res;
}

double Solver::energy(const std::vector<double>& params, const std::vector<double>& transform)
{
	double sum_ = 0;
//	#pragma omp parallel(8)
//	#pragma omp for
	for (size_t i = 0; i < lattice.atoms.size(); ++i)
	{
		double sum_b = 0, sum_r = 0;
		for (size_t j = 0; j < lattice.atoms.size(); ++j)
		{
			if (i == j)
				continue;
			double distance = lattice.distance(i, j, transform);
			if (distance != 0)
			{
				sum_r += ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1)));
				sum_b += params[5] * params[5] * exp(-2 * params[4] * (distance / params[2] - 1));
			}
		}
//		#pragma omp atomic
		sum_ += sum_r - sqrt(sum_b);
	}
	return sum_;
}

double Solver::derivative(double energy, const std::vector<double>& params, const std::vector<double>& pos, const std::vector<double>& neg)
{
	double positiveEnergy = this->energy(params, pos);
	double negativeEnergy = this->energy(params, neg);

	return (positiveEnergy - 2 * energy + negativeEnergy) / (1e-2 * 1e-2);
}
