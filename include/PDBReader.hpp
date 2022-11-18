#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <cctype>
#include <algorithm>

class PDBReader {
    public:
        PDBReader ();
        ~PDBReader ();

        void read_file (const std::string& filename);
        std::vector<std::string> get_atom_names () const;
        std::vector<std::string> get_elements () const;
        std::vector<std::string> get_residue_names () const;
        std::vector<int> get_chain_ids () const;
        std::vector<int> get_residue_numbers () const;
        std::vector<std::array<float,3>> get_atom_coordinates () const;
        std::vector<float> get_atom_temp_factors () const;
    private:
        //changes phosphate names from O1P/O2P to OP1/OP2
        std::string del_ws (const std::string& s) {
            std::string out = s;
            //out.erase(std::remove_if(out.begin(), out.end(),
            //            [] (char c) {return static_cast<unsigned char>(c == '\r'||c=='\t'||c==' '||c=='\n');}), out.end());
            out.erase(std::remove_if(out.begin(), out.end(), ::isspace), out.end());
            return out;
        }
        void change_OP_names ();
        std::vector<int> atom_serial_numbers;
        std::vector<std::string> atom_names;
        std::vector<char> atom_alt_locations;
        std::vector<std::string> atom_residue_names;
        std::vector<int> atom_chain_ids;
        std::vector<int> residue_sequence_numbers;
        std::vector<char> insert_residue_codes;
        std::vector<std::array<float,3>> atom_coords;
        std::vector<float> occupancy;
        std::vector<float> atom_temp_factors;
        std::vector<std::string> atom_elements;
        std::vector<std::string> atom_charges;
};

