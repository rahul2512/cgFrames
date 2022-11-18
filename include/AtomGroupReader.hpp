#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>

#include "atom_group.hpp"

class AtomGroupReader {
    public:
        AtomGroupReader () {}

        AtomGroupReader (const std::string& filename) {
            read_file(filename);
        }

        ~AtomGroupReader () {}

        void read_file (const std::string& filename);

        std::vector<atom_group> get_atom_groups ();

    private:

        std::vector<atom_group> atom_groups;
};
