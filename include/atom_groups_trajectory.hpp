#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <array>

struct atom_groups_trajectory {
    std::vector<std::string> atom_group_names;
    std::vector<std::vector<std::string>> atom_names;
    std::vector<std::vector<int>> atom_indices; //indices of in original data
    //|frames   |atom group |atoms      |coordinates
    std::vector<std::vector<std::vector<std::array<float,3>>>> coordinates;
    std::vector<int> chain_ids;
    std::vector<int> residue_numbers;
    std::vector<std::string> residues;
};
