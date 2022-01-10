#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Atom.h"
#include <set>

struct Lattice
{
	std::vector<Atom> atoms;
	Point period;
	double multiplier;
public:
	Lattice(std::string file_name, Point period, double m);
	Lattice(const Lattice& lat);
	Lattice() { multiplier = 4.085; };
	double dist(size_t i, size_t j,const std::vector<double>& transform , bool cut = false);
	Atom& operator [](int i);
	static Lattice GenerateLattice(Atom::ElementType type, double multiplier);
	Lattice(const std::vector<Atom>& atoms, Point period, double multiplier);
};

