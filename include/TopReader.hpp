#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

class TopReader {
    public:
        TopReader ();
        ~TopReader ();

        void read_file (const std::string& filename);
        //================================================================================
        //  Getters
        //================================================================================
        std::vector<std::string> get_atom_names () const;
        std::vector<std::string> get_residue_names () const;
        std::vector<int> get_residue_numbers () const;
        std::vector<std::string> get_elements () const;
        std::vector<int> get_chain_ids () const;

        void dispose_of_hydrogen();
        void dispose_of_ions();
    private:
        void change_OP_names ();
        std::pair<int,int> determine_format(const std::string& format_string) const;

        //================================================================================
        //  Reading members:
        //================================================================================
        int read_number_of_atoms (const std::string& filename);
        int read_number_of_residues (const std::string& filename);
        std::vector<std::string> read_atom_names (const std::string& filename);
        std::vector<std::string> read_residue_names (const std::string& filename);
        std::vector<int> determine_residue_numbers (const std::string& filename);
        std::vector<std::string> read_element_names (const std::string& filename);
        std::vector<int> determine_chain_IDs (const std::string& filename);

        //================================================================================
        //  Members:
        //================================================================================

        std::vector<std::string> atom_names;
        std::vector<std::string> atom_residue_names;
        std::vector<int> atom_residue_numbers;
        std::vector<std::string> elements;
        std::vector<int> chain_id;
        int number_of_atoms = 0;
        int number_of_residues = 0;
};
