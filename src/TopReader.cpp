#include "TopReader.hpp"

TopReader::TopReader () {}
TopReader::~TopReader () {}

//================================================================================
//  Getters
//================================================================================

std::vector<std::string> TopReader::get_atom_names () const
{
    return atom_names;
}

std::vector<std::string> TopReader::get_residue_names () const
{
    return atom_residue_names;
}

std::vector<int> TopReader::get_residue_numbers () const
{
    return atom_residue_numbers;
}

std::vector<std::string> TopReader::get_elements () const
{
    return elements;
}

std::vector<int> TopReader::get_chain_ids () const
{
    return chain_id;
}

//================================================================================
//  Reading file
//================================================================================

void TopReader::read_file (const std::string& filename)
{
    std::cout << "Reading topology file" << std::endl;

    number_of_atoms = read_number_of_atoms(filename);
    number_of_residues = read_number_of_residues(filename);
    std::cout << "\treading atom names...";
    atom_names = read_atom_names(filename);
    std::cout << "done" << std::endl;

    std::cout << "\treading residue names...";
    atom_residue_names = read_residue_names(filename);
    std::cout << "done" << std::endl;

    std::cout << "\treading residue numbers...";
    atom_residue_numbers = determine_residue_numbers(filename);
    std::cout << "done" << std::endl;

    std::cout << "\treading elements...";
    elements = read_element_names (filename);
    std::cout << "done" << std::endl;

    std::cout << "\tdetermining chain IDs...";
    chain_id = determine_chain_IDs (filename);
    std::cout << "done" << std::endl;


    //Can add other things, like atom type, atomic number, mass, etc
    change_OP_names();

    /*
    std::cout << "atom_names size" << std::endl;
    std::cout << atom_names.size() << std::endl;
    std::cout << "atom_residue_names size" << std::endl;
    std::cout << atom_residue_names.size() << std::endl;
    std::cout << "atom_residue_numbers size" << std::endl;
    std::cout << atom_residue_numbers.size() << std::endl;
    std::cout << "elements size" << std::endl;
    std::cout << elements.size() << std::endl;
    std::cout << "chains size" << std::endl;
    std::cout << chain_id.size() << std::endl;
    */
}

int TopReader::read_number_of_atoms (const std::string& filename)
{
    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }
    std::string line = "";
    std::stringstream ss;

    //Get number of atoms:
    while (std::getline(top_file, line)) {
        if (line.find("FLAG POINTERS") != std::string::npos) {
            std::getline(top_file, line);
            break;
        }
    }
    std::getline(top_file, line); //Read line with number of atoms
    ss << line;
    ss >> number_of_atoms;
    return number_of_atoms;
}

int TopReader::read_number_of_residues (const std::string& filename)
{
    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }
    top_file.seekg(0, top_file.beg); //return to top of file
    std::string line = "";
    std::stringstream ss;

    //Get number of atoms:
    while (std::getline(top_file, line)) {
        if (line.find("FLAG POINTERS") != std::string::npos) {
            std::getline(top_file, line);
            std::getline(top_file, line);
            break;
        }
    }
    std::getline(top_file, line); //Read line with number of residues in it
    ss << line;
    ss >> number_of_residues;
    ss >> number_of_residues;
    return number_of_residues;
}

std::vector<std::string> TopReader::read_atom_names (const std::string& filename)
{
    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }
    //Get atom indices and names:
    std::vector<std::string> atoms;
    std::string line;
    while (std::getline(top_file, line)) {
        if (line.find("FLAG ATOM_NAME") != std::string::npos) {
            break;
        }
    }
    
    //Get format
    while (std::getline(top_file, line)) {
        if (line.find("FORMAT") != std::string::npos) {
            break;
        }
    }

    auto format = determine_format(line);
    int atoms_per_line = format.first;
    int chars_per_atom = format.second;

    //Get the actual atom names:
    int indexed_atoms = 0;
    while (indexed_atoms < number_of_atoms) {
        std::getline(top_file, line);
        std::string name = "";
        for (int a = 0; a < atoms_per_line && indexed_atoms < number_of_atoms; a++) {
            std::string name = line.substr(a*chars_per_atom,chars_per_atom);
            if (name.find(' ') != std::string::npos)
                name.erase(name.find(' '));
            atoms.push_back(name);
            indexed_atoms++;
        }
    }
    return atoms;
}    

std::vector<std::string> TopReader::read_residue_names (const std::string& filename)
{
    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }
    std::string line = "";
    std::stringstream ss;

    //Get atom indices and names:
    std::vector<std::string> residues;
    while (std::getline(top_file, line)) {
        if (line.find("FLAG RESIDUE_LABEL") != std::string::npos) {
            break;
        }
    }
    
    //Get format
    while (std::getline(top_file, line)) {
        if (line.find("FORMAT") != std::string::npos) {
            break;
        }
    }

    auto format = determine_format(line);
    int res_per_line = format.first;
    int chars_per_res = format.second;
    //std::cout << "format:" << std::endl;
    //std::cout << res_per_line << '\t' << chars_per_res << std::endl;

    //Get the actual residue names:
    int indexed_res = 0;
    while (indexed_res < number_of_residues) {
        std::getline(top_file, line);
        std::string name = "";
        for (int a = 0; a < res_per_line && indexed_res < number_of_residues; a++) {
            std::string name = line.substr(a*chars_per_res,chars_per_res);
            if (name.find(' ') != std::string::npos)
                name.erase(name.find(' '));
            residues.push_back(name);
            indexed_res++;
        }
    }

    std::vector<int> residue_indices = determine_residue_numbers(filename);
    
    //couple both
    std::vector<std::string> residue_labels(number_of_atoms);
    //std::cout<<residue_labels.size() << std::endl;
    for (int i = 0; i < number_of_atoms; i++) {
        residue_labels.at(i) = residues[residue_indices[i]];
        //residue_labels.at(residue_indices[i]) = residues[i];
        //residue_labels.at(residue_indices[i]-1) = residues[i];
    }
    //now fill in the rest:
    //std::string current_residue = "";
    //for (int i = 0; i < number_of_atoms; i++) {
    //    if (residue_labels[i] != "")
    //        current_residue = residue_labels[i];
    //    if (residue_labels[i] == "")
    //        residue_labels[i] = current_residue;
    //}


    //for (int i = 0; i < number_of_atoms; i++) 
    //    std::cout << i << '\t' << residue_labels.at(i) << std::endl;

    return residue_labels;
}

std::vector<int> TopReader::determine_residue_numbers (const std::string& filename)
{
    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }

    std::string line = "";
    std::stringstream ss;

    //Find the residue pointers:
    std::vector<int> residue_indices;
    while (std::getline(top_file, line)) {
        if (line.find("FLAG RESIDUE_POINTER") != std::string::npos) {
            break;
        }
    }
    
    //Get format
    while (std::getline(top_file, line)) {
        if (line.find("FORMAT") != std::string::npos) {
            break;
        }
    }

    auto format = determine_format(line);
    int res_per_line = format.first;
    int chars_per_res = format.second;
    //std::cout << "format: " << format.first << '\t' << format.second << std::endl;

    //Get the actual residue indices:
    int indexed_res = 0;
    while (indexed_res < number_of_residues) {
        std::getline(top_file, line);
        //std::cout << "line: " << line << std::endl;
        for (int a = 0; a < res_per_line && indexed_res < number_of_residues; a++) {
            std::string index = line.substr(a*chars_per_res,chars_per_res);
            //index.erase(std::remove_if( index.begin(), index.end(), 
            //    [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}), index.end() );
            index.erase(std::remove_if(index.begin(), index.end(), ::isspace), index.end());


            //std::cout << "index: " << index << std::endl;
            residue_indices.push_back(std::stoi(index)-1);
            indexed_res++;
        }
    }
    std::vector<int> residue_numbers(number_of_atoms);

    int id = 1;
    for (int i = 0; i < number_of_residues; i++) {
        residue_numbers.at(residue_indices[i]) = id++;
    }

    //now fill in the rest:
    int current_residue = 0;
    for (int i = 0; i < number_of_atoms; i++) {
        if (residue_numbers[i] != 0)
            current_residue = residue_numbers[i];
        if (residue_numbers[i] == 0)
            residue_numbers[i] = current_residue;
    }

    for (int i = 0; i < number_of_atoms; ++i)
        residue_numbers[i] -= 1;

    //for (auto& r : residue_numbers)
    //    std::cout << r << std::endl;

    return residue_numbers;
}

std::vector<std::string> TopReader::read_element_names (const std::string& filename)
{

    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }
    std::string line = "";
    std::stringstream ss;

    //Get atom indices and names:
    std::vector<std::string> elements;
    while (std::getline(top_file, line)) {
        if (line.find("FLAG AMBER_ATOM_TYPE") != std::string::npos) {
            break;
        }
    }
    
    //Get format
    while (std::getline(top_file, line)) {
        if (line.find("FORMAT") != std::string::npos) {
            break;
        }
    }

    auto format = determine_format(line);
    int atoms_per_line = format.first;
    int chars_per_atom = format.second;

    //Get the actual element names:
    int indexed_atoms = 0;
    while (indexed_atoms < number_of_atoms) {
        std::getline(top_file, line);
        std::string name = "";
        for (int a = 0; a < atoms_per_line && indexed_atoms < number_of_atoms; a++) {
            std::string name = line.substr(a*chars_per_atom,chars_per_atom);
            if (name.find(' ') != std::string::npos)
                name.erase(name.find(' '));
            elements.push_back(name);
            indexed_atoms++;
        }
    }
    return elements;
}

std::vector<int> TopReader::determine_chain_IDs (const std::string& filename)
{
    std::ifstream top_file(filename, std::ios_base::in);

    if (!top_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }

    std::string line = "";
    std::stringstream ss;

    //Get atom indices and names:
    auto chain_IDs = std::vector<int>(number_of_atoms);
    while (std::getline(top_file, line)) {
        if (line.find("FLAG ATOMS_PER_MOLECULE") != std::string::npos) {
            break;
        }
    }
    
    //Get format
    while (std::getline(top_file, line)) {
        if (line.find("FORMAT") != std::string::npos) {
            break;
        }
    }

    auto format = determine_format(line);
    int numbers_per_line = format.first;
    int chars_per_number = format.second;

    //read it:
    int indexed_atoms = 0;
    int id = 0;
    while (indexed_atoms < number_of_atoms) {
        std::getline(top_file, line);
        std::string name = "";
        for (int a = 0; a < numbers_per_line && indexed_atoms < number_of_atoms; a++) {
            std::string index = line.substr(a*chars_per_number,chars_per_number);
            index.erase(std::remove_if(index.begin(), index.end(), ::isspace), index.end());
            int amount = std::stoi(index);
            for (int i = indexed_atoms; i < indexed_atoms+amount; ++i)
                chain_IDs.at(i) = id;
            indexed_atoms += amount;
            id++;
        }
    }

    return chain_IDs;
}

//returns pair with {atoms_per_line, chars_per_atom}
std::pair<int,int> TopReader::determine_format(const std::string& format_string) const
{
    int atoms_per_line = 0;
    int chars_per_atom = 0;
    std::string format = "";
    bool in_paren = false;
    std::string line;

    for (int c = 0; c < format_string.size(); c++) {
        if (format_string.at(c) == '(') {
            in_paren = true;
            continue; //skip to next c
        }
        if (format_string.at(c) == ')') {
            in_paren = false;
            break;
        }
        if (in_paren)
            format.push_back(format_string[c]);
    }
    std::string atoms_format = format;
    std::string chars_format = format;
    int format_index = -1;
    for (int i = 0; i < format.size(); i++) {
        if (isalpha(format.at(i))) {
            format_index = i; 
            break;
        }
    }
    atoms_per_line = std::stoi(atoms_format.erase(format_index));
    chars_per_atom = std::stoi(chars_format.substr(format_index+1, format.size()-format_index));

    return {atoms_per_line, chars_per_atom};
}


void TopReader::dispose_of_hydrogen ()
{
    for (int i = atom_names.size()-1; i >= 0; i--) {
        if (atom_names[i][0] == 'H') {
            atom_names.erase(atom_names.begin() + i);
        }
    }
}

void TopReader::dispose_of_ions ()
{
    for (int i = atom_names.size()-1; i >= 0; i--) {
        if (atom_names[i] == "K+" || atom_names[i] == "Na+" || atom_names[i] == "Cl-" ) {
            atom_names.erase(atom_names.begin() + i);
        }
    }
}

void TopReader::change_OP_names ()
{
    for (int i = 0; i < atom_names.size(); ++i) {
//	std::cout << atom_names.size() << "atom_names.size()" <<std::flush;
        if (atom_names[i] == "O1P")
            atom_names[i] = "OP1";
        else if (atom_names[i] == "O2P")
            atom_names[i] = "OP2";
//        else if (atom_names[1] == "O5'")
//            atom_names[1] = "5O'";
//	else if (atom_names[764] == "O5'")
//            atom_names[764] = "5O'";

    }
}

