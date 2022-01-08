#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Atom.h"

struct Lattice
{
	std::vector<Atom> atoms;
	Point period;
	double multiplier;
public:
	Lattice(std::string file_name, Point transform, double m);
	Lattice(const Lattice& lat);
	Lattice() { multiplier = 1; };
	double dist(size_t i, size_t j,const std::vector<double>& period, bool cut);
	Atom& operator [](int i);
};

