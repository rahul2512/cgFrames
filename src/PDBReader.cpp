#include "PDBReader.hpp"

PDBReader::PDBReader () {}
PDBReader::~PDBReader () {}

void PDBReader::read_file (const std::string& filename)
{
    std::ifstream pdb_file(filename, std::ios_base::in);
    if (!pdb_file.is_open()) {
        std::cout << "topology file not opened" << std::endl;
        exit(1);
    }
    std::string line = "";

    //skip all non-atom parts of the PDB
    while (std::getline(pdb_file, line)) {
        if (line.substr(0,6) == "ATOM  ")
            break;
    }

    do {
        //Check if we reached the end:
        if (line.substr(0,3) == "END")
            break;
        if (line.substr(0,3) == "REM")
            continue;
        if (line.substr(0,3) == "HEA")
            continue;
        if (line.substr(0,3) == "TIT")
            continue;	
        if (line.substr(0,3) == "COM")
            continue;
        if (line.substr(0,3) == "AUT")
            continue;
        if (line.substr(0,3) == "JRN")
            continue;
        if (line.substr(0,3) == "REV")
            continue;
        if (line.substr(0,3) == "EXP")
            continue;
        if (line.substr(0,3) == "KEY")
            continue;
        if (line.substr(0,3) == "DBR")
            continue;
        if (line.substr(0,3) == "SEQ")
            continue;
        if (line.substr(0,3) == "FOR")
            continue;
        if (line.substr(0,3) == "ORI")
            continue;
        if (line.substr(0,3) == "CRY")
            continue;
        if (line.substr(0,3) == "SCA")
            continue;
        if (line.substr(0,3) == "HET")
            continue;
	if (line.substr(0,3) == "TER")
           continue;
        if (line.substr(0,3) == "MAS")
           continue;
	
        //If not, parse the atom data

        //PDB lists the ATOM records as follows:
        //col       description
        //===========================
        //1-6       "ATOM  "
        //7-11      atom serial number
        //13-16     atom name
        //17        Alt location ID
        //18-20     residue name
        //22        chain ID
        //23-26     Residue sequence number
        //27        Code for insertion of residues
        //31-38     x-coord of atom
        //39-46     y-coord of atom
        //47-54     z-coord of atom
        //55-60     occupancy
        //61-66     temperature factor
        //77-78     element symbol
        //79-80     atom charge
//        std::cout << std::endl;
//        std::cout << "0         0         0         0         0         0         0" << std::endl;
//        std::cout << "012345678901234567890123456789012345678901234567890123456789" << std::endl;
//        std::cout << line << std::endl;
//        std::cout << "line size: " << line.size() << std::endl;

//        std::cout << "atom number:" << '\t';
//        std::cout << del_ws(line.substr(6,5)) << std::endl;

        atom_serial_numbers.push_back(std::stoi(del_ws(line.substr(6,5))));
//        std::cout << "atom name:" << '\t';
//        std::cout << del_ws(line.substr(12,4)) << std::endl;

        atom_names.push_back(del_ws(line.substr(12,4)));

//        std::cout << "atom alt loc:" << '\t';
//        std::cout << line.substr(16,1) << std::endl;

        atom_alt_locations.push_back(line.at(16));

//        std::cout << "atom residue:" << '\t';
//        std::cout << del_ws(line.substr(17,3)) << std::endl;

        atom_residue_names.push_back(del_ws(line.substr(17,3)));

//        std::cout << "atom chain id:" << '\t';
//        std::cout << line.substr(21,1) << std::endl;

        atom_chain_ids.push_back(line.at(21));   //rahul

//        std::cout << "residue seq no:" << '\t';
//        std::cout << del_ws(line.substr(23,3)) << std::endl;

        residue_sequence_numbers.push_back(std::stoi(del_ws(line.substr(23,3))));

//        std::cout << "inert residue code:" << '\t';
//        std::cout << del_ws(line.substr(27,1)) << std::endl;

        insert_residue_codes.push_back(line.at(27));

//        std::cout << "atom coords:" << '\t';
//        std::cout << del_ws(line.substr(30,8)) << '\t';
//        std::cout << del_ws(line.substr(38,8)) << '\t';
//       std::cout << del_ws(line.substr(46,8)) << std::endl;

        atom_coords.push_back({std::stof(del_ws(line.substr(30,8))),
                               std::stof(del_ws(line.substr(38,8))),
                               std::stof(del_ws(line.substr(46,8)))});

//        std::cout << "occupancy:" << '\t';
//        std::cout << del_ws(line.substr(54,5)) << std::endl;

        occupancy.push_back(std::stof(del_ws(line.substr(54,5))));

//        std::cout << "temp factor:" << '\t';
//        std::cout << del_ws(line.substr(60,5)) << std::endl;

        atom_temp_factors.push_back(std::stof(del_ws(line.substr(60,5))));

//        std::cout << "atom element:" << '\t';
//        std::cout << del_ws(line.substr(76,2)) << std::endl;

        atom_elements.push_back(del_ws(line.substr(76,2)));

//        std::cout << "atom charge:" << '\t';

        if (line.size() > 78)
            atom_charges.push_back(del_ws(line.substr(78,2)));
        else atom_charges.push_back("");
    } while (std::getline(pdb_file, line));

    change_OP_names();

    pdb_file.close();
}

std::vector<std::string> PDBReader::get_atom_names () const
{
    return atom_names;
}

std::vector<std::string> PDBReader::get_elements () const
{
    return atom_elements;
}

std::vector<std::string> PDBReader::get_residue_names () const
{
    return atom_residue_names;
}

std::vector<int> PDBReader::get_chain_ids () const
{
    return atom_chain_ids;
}

std::vector<int> PDBReader::get_residue_numbers () const
{
    return residue_sequence_numbers;
}

std::vector<std::array<float,3>> PDBReader::get_atom_coordinates () const
{
    return atom_coords;
}

std::vector<float> PDBReader::get_atom_temp_factors () const
{
    return atom_temp_factors;
}

void PDBReader::change_OP_names ()
{
    for (int i = 0; i < atom_names.size(); ++i) {
        if (atom_names[i] == "O1P")
            atom_names[i] = "OP1";
        else if (atom_names[i] == "O2P")
            atom_names[i] = "OP2";
    }
}
