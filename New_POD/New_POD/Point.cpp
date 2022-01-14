#include "Point.h"

Point::Point(double x, double y, double z): x(x), y(y), z(z)
{}

Point::Point(const std::vector<double>& vec)
{
	if (vec.size() == 3)
		x = vec[0], y = vec[1], z = vec[2];
	else
		x = y = z = 0;
}

double& Point::operator[](size_t i)
{
	switch (i)
	{
	case 0: return x;
	case 1: return y;
	case 2: return z;
	}
}

double Point::norm(const Point& p)
{
	return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

Point::Point()
{
	x = y = z = 0;
}