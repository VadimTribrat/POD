#include <omp.h>
#include "Solver.h"

Solver::Solver(const Lattice& lat, double e):lattice(lat), E_cohA(e), E_cohB(0)
{
	fittedAB = fittedBB = fittedAA = std::vector<double>();
}


std::vector<double> Solver::calculateBB(const std::vector<double>& params)
{
	std::vector<double> res;
	lattice.serMultiplier(calculateLatticeConstant(params, 3.95, 4.2, 5, 5));
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
			if (i != j)
			{
				double distance = lattice.distance(i, j, transform);
				if (distance != 0)
				{
					if (lattice.atoms[i].first == Lattice::ElementType::B && lattice.atoms[i].first == lattice.atoms[j].first)
					{
						if (fittedBB.empty())
						{
							sum_r += ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1)));
							sum_b += params[5] * params[5] * exp(-2 * params[4] * (distance / params[2] - 1));
						}
						else
						{
							sum_r += ((fittedBB[0] * (distance - fittedBB[2]) + fittedBB[1]) * exp(-fittedBB[3] * (distance / fittedBB[2] - 1)));
							sum_b += fittedBB[5] * fittedBB[5] * exp(-2 * fittedBB[4] * (distance / fittedBB[2] - 1));
						}
					}
					else if (lattice.atoms[i].first != lattice.atoms[j].first)
					{
						if (fittedAB.empty())
						{
							sum_r += ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1)));
							sum_b += params[5] * params[5] * exp(-2 * params[4] * (distance / params[2] - 1));
						}
						else
						{
							sum_r += ((fittedAB[0] * (distance - fittedAB[2]) + fittedAB[1]) * exp(-fittedAB[3] * (distance / fittedAB[2] - 1)));
							sum_b += fittedAB[5] * fittedAB[5] * exp(-2 * fittedAB[4] * (distance / fittedAB[2] - 1));
						}
					}
					else
					{
						if (fittedAA.empty())
						{
							sum_r += ((params[0] * (distance - params[2]) + params[1]) * exp(-params[3] * (distance / params[2] - 1)));
							sum_b += params[5] * params[5] * exp(-2 * params[4] * (distance / params[2] - 1));
						}
						else
						{
							sum_r += ((fittedAA[0] * (distance - fittedAA[2]) + fittedAA[1]) * exp(-fittedAA[3] * (distance / fittedAA[2] - 1)));
							sum_b += fittedAA[5] * fittedAA[5] * exp(-2 * fittedAA[4] * (distance / fittedAA[2] - 1));
						}
					}
				}
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

double Solver::calculateLatticeConstant(const std::vector<double>& params, double a_min, double a_max, int n, int itersCount) {
	std::pair<double, double> min(0, 1e30);
	double width;

	int iter = 0;
	do {
		width = a_max - a_min;
		min = std::make_pair(0, 1e8);

		for (int i = 0; i < n; ++i) {
			double point = a_min + i * width / n;

			lattice.serMultiplier(point);
			double value = energy(params);

			if (value < min.second) {
				min = std::make_pair(point, value);
			}
		}

		width /= n;
		a_min = min.first - width / 2;
		a_max = min.first + width / 2;

		++iter;
	} while (iter < itersCount);

	return min.first;
}

std::vector<double> Solver::calculateAB(const std::vector<double>& params)
{
	double E_B = energy(params);
	lattice[0].first = Lattice::ElementType::A;
	auto E_AB = energy(params);
	lattice[0].first = Lattice::ElementType::B;
	return std::vector<double>{E_AB - E_B - E_cohA + E_cohB};
}

std::vector<double> Solver::calculateAA(const std::vector<double>& params)
{
	auto lat = Lattice(Point(3, 3, 6), 4.085, true);
	double E_surf = energy(params);
	//lat.print();
	lat[12].first = Lattice::ElementType::A;
	double E_adaatom = energy(params);
	lat[68].first = Lattice::ElementType::A;
	double E_dim = energy(params);
	lat[12].first = lat[68].first = Lattice::ElementType::B;
	double E_indim = E_dim - E_surf - 2 * (E_adaatom - E_surf);
	lat.atoms.push_back(std::make_pair(Lattice::ElementType::A, Point(0.5, 0, 3.5)));
	E_adaatom = energy(params);
	lat.atoms.push_back(std::make_pair(Lattice::ElementType::A, Point(0, 0.5, 3.5)));
	
	E_dim = energy(params);
	double E_ondim = E_dim - E_surf - 2 * (E_adaatom - E_surf);
	return std::vector<double>{E_indim, E_ondim};
}