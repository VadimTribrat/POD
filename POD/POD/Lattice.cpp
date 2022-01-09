#include "Lattice.h"

Lattice::Lattice(std::string file_name, Point period, double m) : period(period), multiplier(m)
{
    std::ifstream in(file_name);
    std::string line;
    std::set<std::vector<double>> set_;
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
    for (size_t i = 0; i<atoms.size(); ++i)
    {
        std::cout << i << " " << atoms[i] << "\n";
    }
    std::cout << "\n";
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

    if (transform.size() == 3) {
        dir = Point::multiply(dir, Point(transform));
    }
    else if (transform.size() == 5) {
        dir = Point(dir.x * transform[0] + dir.y * transform[1],
            dir.x * transform[2] + dir.y * transform[3],
            dir.z * transform[4]);
    }

    double norm = Point::norm(dir);
    if (norm >= sqrt(3) && true) {
        return 0;
    }

    return norm*multiplier;
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