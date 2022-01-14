#pragma once
#include <iostream>
#include <vector>

struct Point
{
	double x, y, z;
	Point(double x, double y, double z);
	Point(const std::vector<double>& vec);
	Point(const Point& p) = default;
	Point& operator=(const Point& p) = default;
	Point();
	double& operator[](size_t i);
	static double norm(const Point& p);
	friend Point operator*(const Point& p1, const Point& p2)
	{
		return Point(p1.x * p2.x, p1.y * p2.y, p1.z * p2.z);
	}
	friend Point operator+(const Point& p1, const Point& p2)
	{
		return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
	}
	friend Point operator-(const Point& p1, const Point& p2)
	{
		return Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
	}
	friend std::ostream& operator<< (std::ostream& out, const Point& p)
	{
		out << p.x << " " << p.y << " " << p.z;
		return out;
	}
};