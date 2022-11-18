#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>

#include <netcdf>

#include "atom_group.hpp"

class TrajReader {
    public:
        TrajReader () {}
        ~TrajReader () {}

        void read_file (const std::string& filename);
        void dispose_of_hydrogen (const std::vector<std::string>& atom_names);
        void dispose_of_ions (const std::vector<std::string>& atom_names);
        std::vector<std::vector<std::array<float,3>>> get_raw_trajectories ();

        //atom groups, frames,  atoms,       coordinates
        //std::vector<atom_group_traj> get_atom_group_trajectories (const std::vector<std::string>& atom_names, const std::vector<atom_group>& ideal_bases);
    private:
        std::vector<std::vector<std::array<float,3>>> parse_raw_traj (int n_atoms, int n_dims, int n_frames, float* data);
        int amount_frames = 0;
        int amount_atoms = 0;
        int amount_dims = 0;
        std::vector<std::vector<std::array<float,3>>> raw_trajectories;
};
