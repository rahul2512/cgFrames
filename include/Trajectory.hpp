#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <iterator>
#include <algorithm>

#include "Topology.hpp"
#include "PDBReader.hpp"
#include "TopReader.hpp"
#include "TrajReader.hpp"

enum traj_query {
    DNA_bases,
    DNA_backbone,
    ions,
    water,
    DNA_hydrogen,
    DNA_phosphates,
    DNA_sugar
};

class Trajectory {
    public:
        //==============================================================================
        //  Con/destructors
        //==============================================================================

        Trajectory () {}
        ~Trajectory () {}

        //==============================================================================
        //Load data from file:
        //==============================================================================

        /**
         * Load trajectory and topology file. If trajectory file is .nc, topology file is necessary.
         * If trajectory file is pdb, topology file can be left blank to use the same file for the
         * topology, or if one is specified, use that file.
         * Topology files can be .top or .pdb.
         */
        void load (const std::string& trajectory_filename, const std::string& topology_filename= "");

        //==============================================================================
        //  Queries/ Getters
        //==============================================================================

        int get_number_of_md_frames () {
            return n_frames;
        }

        std::vector<int> select (const traj_query& q) const;
        std::vector<int> select_inverse (const std::vector<int>& indices) const;
        std::vector<int> select_element (const std::string& element) const;
        std::vector<int> select_atom_with_name (const std::string& atom_name) const;

        void remove (const std::vector<int>& indices = std::vector<int>());

        std::vector<std::string> get_names (const std::vector<int>& indices = std::vector<int>()) const;
        std::vector<std::string> get_residues (const std::vector<int>& indices = std::vector<int>()) const;
        std::vector<int> get_residue_numbers (const std::vector<int>& indices = std::vector<int>()) const;
        std::vector<std::string> get_elements (const std::vector<int>& indices = std::vector<int>()) const;
        std::vector<int> get_chain_ids (const std::vector<int>& indices = std::vector<int>()) const;

        std::vector<std::vector<float*>> get_coordinates (const std::vector<int>& indices = std::vector<int>()) const;
        std::vector<std::vector<std::array<float,3>>> get_coordinates_stl (const std::vector<int>& indices = std::vector<int>()) const;
    private:

        //==============================================================================
        //  Read different data formats
        //==============================================================================

        Topology read_PDB_topology (const std::string& topology_filename) const;
        Topology read_TOP_topology (const std::string& topology_filename) const;
        std::vector<std::vector<std::array<float,3>>> read_netCDF_trajectory (const std::string& trajectory_filename) const;
        std::vector<std::vector<std::array<float,3>>> read_PDB_trajectory (const std::string& trajectory_filename) const;
        std::vector<std::vector<std::array<float,3>>> parse_raw_pdbtraj (int n_atoms, int n_dims, int n_frames, std::vector<std::vector<std::array<float,3>>> data);

        //==============================================================================
        //  Member data
        //==============================================================================

        int n_frames;
        int n_atoms;
        //frames         atoms  coords
        std::vector<std::vector<std::array<float,3>>> coordinates;
	std::vector<std::vector<std::array<float,3>>> coor;
        //std::vector<std::vector<float*>> coordinates;
        Topology topology;
};
