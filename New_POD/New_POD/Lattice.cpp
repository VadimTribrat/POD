#include "Lattice.h"

Lattice::Lattice(const std::string& file_name, const Point& period, double m) : period(period), multiplier(m)
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
                    atoms.push_back(std::make_pair(Lattice::ElementType::A,
                        Point(std::stod(temp[1]), std::stod(temp[2]), std::stod(temp[3]))));
                else
                    atoms.emplace_back(std::make_pair(Lattice::ElementType::B,
                        Point(std::stod(temp[1]), std::stod(temp[2]), std::stod(temp[3]))));
        }
        in.close();
    }
    else
        throw std::exception("");
}

Lattice::Lattice(int d, double m)
{
    period = Point(d, d, d);
    multiplier = m;
    std::vector<Point> basis = {
        Point(0, 0, 0),
        Point(0.5, 0, 0.5),
        Point(0.5, 0.5, 0),
        Point(0, 0.5, 0.5)
    };

    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            for (int k = 0; k < d; ++k) {
                for (auto& val : basis) {
                    atoms.push_back(std::make_pair(Lattice::ElementType::A, val + Point(i, j, k)));
                }
            }
        }
    }
}

void Lattice::print()
{
    for (size_t i = 0; i < atoms.size(); ++i)
    {
        std::cout << i << " " << atoms[i].second << "\n";
    }
    std::cout << "\n";
}

std::pair<Lattice::ElementType, Point>& Lattice::operator[](int i)
{
    return atoms[i];
}

void Lattice::serMultiplier(double m)
{
    multiplier = m;
}

void Lattice::setPeriod(const Point& p)
{
    period = p;
}

double Lattice::distance(size_t i, size_t j, const std::vector<double>& transform)
{
    auto atom1 = atoms[i].second, atom2 = atoms[j].second;
    auto dis = atom2 - atom1;

    for (size_t i = 0; i < 3; ++i)
        if (dis[i] > period[i] / 2)
            dis[i] -= period[i];
        else if (dis[i] < -period[i] / 2)
            dis[i] += period[i];

    if (transform.size() == 3)
        dis = dis * Point(transform);
    else if (transform.size() == 5)
        dis = Point(dis.x * transform[0] + dis.y * transform[1],
            dis.x * transform[2] + dis.y * transform[3],
            dis.z * transform[4]);

    double norm = Point::norm(dis);

    if (norm > 1.7)
        return 0;
    return norm * multiplier;
}