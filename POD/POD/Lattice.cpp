#include "Lattice.h"

Lattice::Lattice(std::string file_name, Point pos, double m) : period(pos), multiplier(m)
{
    std::ifstream in(file_name);
    std::string line;
    if (in.is_open())
    {
        while (getline(in, line))
        {
            std::vector<std::string> temp;
            std::istringstream input(line);
            while (!input.eof())
            {
                std::string substring;
                input >> substring;
                temp.push_back(substring);
            }
            if (temp.size() == 4)
                if (temp[0] == "Ag")
                    atoms.emplace_back(Atom(Atom::ElementType::Ag,
                        Point(std::stod(temp[1]), std::stod(temp[2]), std::stod(temp[3]))));
                else
                    atoms.emplace_back(Atom(Atom::ElementType::Co,
                        Point(std::stod(temp[1]), std::stod(temp[2]), std::stod(temp[3]))));
        }
    }
    else
        throw std::exception("");
    in.close();
    //for (const auto& val : atoms)
    //    std::cout << val << "\n";
}

double Lattice::dist(size_t i, size_t j, const std::vector<double>& transform, bool cut)
{
    Point dir = atoms[j].pos - atoms[i].pos;
    for (int i = 0; i < 3; ++i) {
        if (dir[i] > period[i] / 2) {
            dir[i] -= period[i];
        }
        else if (dir[i] < -period[i] / 2) {
            dir[i] += period[i];
        }
    }

    double norm = Point::norm(dir);
    //if (norm > 1.7 && cut) {
    //    return std::numeric_limits<double>::max();
    //}

    return norm;
}

Atom& Lattice::operator[](int i)
{
    return atoms[i];
}

Lattice::Lattice(const Lattice& lat)
{
    this->atoms = lat.atoms;
    period = lat.period;
    multiplier = lat.multiplier;
}