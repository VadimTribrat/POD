#pragma once
#include <tuple>
#include <iostream>
#include <Exception>
#include <vector>
#include <cmath>

struct Point
{
	double x, y, z;
public:
	Point(double x, double y, double z);
	Point(const std::vector<double>& vec);
	Point(const Point& p);
	Point();
	double& operator[](int i);
	double operator[](int i) const;
	static double norm(const Point& p);
	static double scalar_product(const Point& p1, const Point& p2);
	static Point multiply(const Point& p1, const Point& p2);
	static bool equal(const Point& p1, const Point& p2);
	friend Point operator-(const Point& p1, const Point& p2)
	{
		return Point(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
	}
	friend Point operator+(const Point& p1, const Point& p2)
	{
		return Point(p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]);
	}
	friend std::ostream& operator<< (std::ostream& out, const Point& p)
	{
		out << p.x << " " << p.y << " " << p.z;
		return out;
	}
};
