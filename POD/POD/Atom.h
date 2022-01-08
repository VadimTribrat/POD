#pragma once
#include "Point.h"

struct Atom
{
public:
	enum class ElementType
	{
		Ag, Co
	};
	enum class ConnectionType
	{
		AA, BB, AB
	};
	static ConnectionType LinkType(const Atom& a1, const Atom& a2);
	Atom(ElementType et, Point p);
	Point pos;
	ElementType type;
	friend std::ostream& operator <<(std::ostream& out, const Atom& atom)
	{
		out << atom.pos;
		return out;
	}
	Atom(const Atom& atom);
};
