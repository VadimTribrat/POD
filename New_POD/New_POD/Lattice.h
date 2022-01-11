#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Point.h"
struct Lattice
{
	enum class ElementType
	{
		A, B
	};
	std::vector<std::pair<ElementType, Point>> atoms;
	Point period;
	double multiplier;
	Lattice(const std::string& file_name, const Point& period, double m);
	Lattice(const Lattice&) = default;
	Lattice(int d, double m);
	Lattice& operator=(const Lattice&) = default;
	std::pair<ElementType, Point>& operator[](int i);
	double distance(size_t i, size_t j, const std::vector<double>& transform = std::vector<double>());
	void setPeriod(const Point& per);
	void serMultiplier(double m);
	void print();
};

