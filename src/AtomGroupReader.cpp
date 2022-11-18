#include "AtomGroupReader.hpp"
#include "utils.hpp"


void AtomGroupReader::read_file (const std::string& filename)
{
    std::ifstream ib_file(filename, std::ios_base::in);
    if (!ib_file.is_open()) {
        std::cout << "ideal bases file not opened" << std::endl;
        exit(1);
    }
    std::string line = "";
    std::stringstream ss;
    std::string name = "";
    std::vector<atom_group> ideal_bases;

    while (std::getline(ib_file, line)) {
        if (line[0] == '#') {
            continue;
        } else {
            //std::cout << line << std::endl;
        //} if (line[0] != '#') {
            auto info = split(line, ' ');
            atom_group group;
            group.name = info.at(0);
            group.n_atoms = std::stoi(info.at(1));
            if (info.size() >= 14) {
                group.R(0,0) = std::stod(info[2]);
                group.R(0,1) = std::stod(info[3]);
                group.R(0,2) = std::stod(info[4]);
                group.R(1,0) = std::stod(info[5]);
                group.R(1,1) = std::stod(info[6]);
                group.R(1,2) = std::stod(info[7]);
                group.R(2,0) = std::stod(info[8]);
                group.R(2,1) = std::stod(info[9]);
                group.R(2,2) = std::stod(info[10]);

                group.r(0) = std::stod(info[11]);
                group.r(1) = std::stod(info[12]);
                group.r(2) = std::stod(info[13]);
            }
            for (int i = 0; i < group.n_atoms; i++) {
                atom a;
                std::getline(ib_file, line);
                auto c = split(line, ' ');
                if (c.size() > 4) {
                    std::cout << "Something wrong with coords in " << group.name << std::endl;
                    std::cout << "Line is: " << line << std::endl;
                    exit(1);
                }
                //std::array<float,3> coord;
                for (int j = 0; j < 3; j++) {
                    a.position[j] = std::stof(c[j]);
                }
                a.name = c.at(3);
                //group.atoms.names.push_back(c.at(3));
                //group.atom_coordinates.push_back(coord);
                group.atoms.push_back(a);
            }
            
            ideal_bases.push_back(group);
        }
    }
    atom_groups = ideal_bases;
    /*
    for (auto& group : atom_groups) {
        std::cout << group.name << "\t" << group.n_atoms << std::endl;

        for (int i = 0; i < group.atom_names.size(); ++i) {
            std::cout << "\t" << group.atom_names[i] << ": " << group.atom_coordinates[i][0] << ' ' << group.atom_coordinates[i][1] << ' ' << group.atom_coordinates[i][2] << ' ' << std::endl;
        }
    }
    */
}

std::vector<atom_group> AtomGroupReader::get_atom_groups ()
{
    return atom_groups;
}
