#include "Point.h"

Point::Point(double xx, double yy, double zz)
{
	x = xx;
	y = yy;
	z = zz;
}

Point::Point(const std::vector<double>& vec)
{
	if (vec.size() == 3)
	{
		x = vec[0], y = vec[1], z = vec[2];
	}
	else
	{
		std::cout << "Point(const std::vector<double>& vec)\n";
		x = y = z = std::numeric_limits<double>::min();
	}
}

double& Point::operator[](int i) 
{
	switch (i)
	{
	case 0: return x; break;
	case 1: return y; break;
	case 2: return z; break;
	}
}

const double& Point::operator[](int i) const
{
	switch (i)
	{
	case 0: return x; break;
	case 1: return y; break;
	case 2: return z; break;
	}
}

double Point::norm(const Point& p)
{
	return std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

double Point::scalar_product(const Point& p1, const Point& p2)
{
	return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

Point Point::multiply(const Point& p1, const Point& p2)
{
	return Point(p1[0] * p2[0], p1[1] * p2[1], p1[2] * p2[2]);
}

bool Point::equal(const Point& p1, const Point& p2)
{
	for (size_t i = 0; i < 3; ++i)
		if (p1[i] != p2[i])
			return false;
	return true;
}

Point::Point(const Point& p)
{

	x = p.x; y = p.y; z = p.z;
}

Point::Point()
{
	x = y = z = 0;
}