#ifndef CGFRAMES_H_
#define CGFRAMES_H_

#include <iostream>
#include <vector>
#include <string>

#include "Trajectory.hpp"
#include "atom_group.hpp"
#include "atom_groups_trajectory.hpp"
#include "AtomGroupReader.hpp"
#include "PDBReader.hpp"
#include "TopReader.hpp"
#include "TrajReader.hpp"
#include "Frame_arma.hpp"
#include "fitter_arma.hpp"
#include "fraWriter.hpp"


class cgFrames {
    public:
        cgFrames () {
            verbose_log = false;
        }
        ~cgFrames () {}

        void read_ideal_bases (const std::string& filename);

        void read_trajectory (const std::string& trajectory_filename,
                              const std::string& topology_filename = "");

        void strip_trajectory ();

        /**
         * Finds the ideal bases in the names of the Topology, then returns
         * the indices to the corresponding atom groups (with vectors to atom
         * indices in Topology) and the names of the ideal bases.
         */
        std::vector<atom_group> find_groups_in_trajectory ();

        /**
         * Creates atom_groups_trajectory from list of atom groups with list
         * of atom indices, atom group names and frame indices.
         */
        std::vector<std::vector<atom_group>> get_atom_group_trajectories (
                                                            const std::vector<int>& frame_indices = std::vector<int>());

        
        void fit_frames_to_trajectory (std::vector<std::vector<atom_group>>& ag_coords);

        //set_fitter (...);
        //set_output (...)

        void set_output_filename (const std::string& filename);

        void set_verbose_log () {
            verbose_log = true;
        }

        FraWriter output;

        //======================================================================
        //  Trajectory getter functions
        //======================================================================
        //...

    private:
        Trajectory trajectory;

        svd_fitter fit;

        //ideal bases:
        std::vector<atom_group> ideal_bases;

        bool verbose_log;
};
#endif//CGFRAMES_H_
