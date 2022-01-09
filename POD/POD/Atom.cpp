#include "Atom.h"

Atom::ConnectionType Atom::LinkType(const Atom& a1, const Atom& a2)
{
	if (a1.type == Atom::ElementType::Ag && a2.type == Atom::ElementType::Ag)
		return Atom::ConnectionType::AA;
	if ((a1.type == Atom::ElementType::Ag || a1.type == Atom::ElementType::Co) && a2.type != a1.type)
		return Atom::ConnectionType::AB;
	return Atom::ConnectionType::BB;
}

Atom::Atom(Atom::ElementType et, Point p): type(et), pos(p)
{}

Atom::Atom(const Atom& atom)
{
	this->pos = atom.pos;
	this->type = atom.type;
}