#pragma once

#include <iostream>
#include <vector>
#include <string>

struct Topology {
    std::vector<std::string> atom_names;
    std::vector<std::string> elements;
    std::vector<std::string> residues;
    std::vector<int> residue_numbers;
    std::vector<int> chain_ids;
};

