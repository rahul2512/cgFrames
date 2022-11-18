#ifndef ATOM_GROUP_H
#define ATOM_GROUP_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <armadillo>

/*
 * IN JSON FORMAT:
 *{
 *  atom_group_name,
 *  atom_group_index,
 *  frame_orientation,
 *  frame_position,
 *  atoms : [
 *      {
 *          atom_name,
 *          atom_position
 *      }, ...
 *  ]
 *}
*/

struct atom {
    std::string name;
    int md_index; //index in topology file from MD
    std::array<float, 3> position;
};

struct atom_group {
    std::string name;
    int index;
    int chain_index;
    int n_atoms;
    std::vector<atom> atoms;
    arma::mat::fixed<3,3> R;
    arma::vec::fixed<3> r;
};

#endif// ATOM_GROUP_H
