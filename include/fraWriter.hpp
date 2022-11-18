#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_map>
//#include "Frame_arma.hpp"
#include "atom_group.hpp"

class FraWriter {
    public:
        FraWriter ();
        ~FraWriter ();
        //configure what ideal base/atom group goes in which file, e.g. bases in .fra, phosphates in .pfra
        //Example  set_output_names ({{"Thymine","out.fra"}, {"Phosphate", "out.pfra"}})
        void set_output_names (const std::vector<std::pair<std::string,std::string>>& names);
        void add_output_name (const std::pair<std::string,std::string>& name);
        void remove_output_name (const std::string& name);
        void set_output_to_dna (const std::string& fra_file_name, const std::string& pfra_file_name);
        void write (std::vector<std::vector<atom_group>>& groups) const;
    private:
        //hash map with atom group names as keys and associated file names for values
        std::unordered_map<std::string,std::string> groups_file_names;
};
