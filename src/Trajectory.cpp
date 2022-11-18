#include "Trajectory.hpp"
#include "utils.hpp"

//==============================================================================
//Load data from file:
//==============================================================================

void Trajectory::load (const std::string& trajectory_filename, const std::string& topology_filename)
{
    //Check if members are empty?
    //
    std::string traj_extension = split(trajectory_filename, '.').back();
    std::string top_extension;
    if (topology_filename != "")
        top_extension = split(topology_filename, '.').back();

    if (traj_extension == "nc") {
        coordinates = read_netCDF_trajectory (trajectory_filename);
        n_frames = coordinates.size();
        n_atoms = coordinates[0].size();
        if (topology_filename == "") {
            std::cout << "netCDF file needs topology file!" << std::endl;
            exit(1);
        } else {
            if (top_extension == "pdb")
                topology = read_PDB_topology(topology_filename);
            else if (top_extension == "top")
                topology = read_TOP_topology(topology_filename);
            else {
                std::cout << "Unsupported topology file format" << std::endl;
                exit(1);
            }
        }
    }
    if (traj_extension == "pdb") {
        coor = read_PDB_trajectory (trajectory_filename);
//        n_frames = 1;
//        n_atoms = coor[0].size();
        if (topology_filename == "") {
            topology = read_PDB_topology(trajectory_filename);
        } else {
            if (top_extension == "pdb")
                topology = read_PDB_topology(topology_filename);
            else if (top_extension == "top"){
                topology = read_TOP_topology(topology_filename);
            } else {
                std::cout << "Unsupported topology file format" << std::endl;
                exit(1);
            	}
            }
//        int n_frames = coor[0].size()/topology.atom_names.size();
        n_frames = coor[0].size()/topology.atom_names.size();
        int n_dims = 3;
        int n_atoms = topology.atom_names.size();
	coordinates = parse_raw_pdbtraj(n_atoms, n_dims, n_frames, coor);
       }
    }

//==============================================================================
//  Queries/ Getters
//==============================================================================

std::vector<int> Trajectory::select (const traj_query& q) const
{
}

//TODO test
std::vector<int> Trajectory::select_inverse (const std::vector<int>& indices) const
{
    std::vector<int> inverse_indices;
    for (int i = 0; i < topology.elements.size(); ++i) {
//	std::cout << "topology.elements.size()" <<topology.elements.size()<<"indices.size     "<<indices.size() <<"\n"<< std::flush ; 
        if (std::find(indices.begin(), indices.end(), i) != indices.end())
            continue;
        else
            inverse_indices.push_back(i);
    }
    return inverse_indices;
}

std::vector<int> Trajectory::select_element (const std::string& element) const
{
    std::vector<int> indices;

    for (int i = 0; i < topology.elements.size(); ++i) {
        //if (topology.elements.at(i) == element)
        if (topology.elements.at(i)[0] == element[0])
            indices.push_back(i);
    }

    return indices;
}

std::vector<int> Trajectory::select_atom_with_name (const std::string& atom_name) const
{
    std::vector<int> indices;

    for (int i = 0; i < topology.atom_names.size(); ++i) {
        //if (topology.elements.at(i) == element)
        if (topology.atom_names.at(i) == atom_name)
            indices.push_back(i);
    }

    return indices;
}

void Trajectory::remove (const std::vector<int>& indices)
{
    //Remove from coordinates
    for (int frame = 0; frame < n_frames; frame++) {
        for (int i = indices.size()-1; i >=0; i--) {
            coordinates[frame].erase(coordinates[frame].begin() + indices[i]);
        }
    }
    //Remove from Topology
    for (int i = indices.size()-1; i >=0; i--) {
        topology.atom_names.erase(topology.atom_names.begin() + indices[i]);
        topology.elements.erase(topology.elements.begin() + indices[i]);
        topology.residues.erase(topology.residues.begin() + indices[i]);
        topology.residue_numbers.erase(topology.residue_numbers.begin() + indices[i]);
        topology.chain_ids.erase(topology.chain_ids.begin() + indices[i]);
    }
}

std::vector<std::string> Trajectory::get_names (const std::vector<int>& indices) const
{
    if (indices.empty())
        return topology.atom_names;

    std::vector<std::string> names;
    //names.reserve(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        names.push_back(topology.atom_names.at(i));
    }
    return names;
}

std::vector<std::string> Trajectory::get_residues (const std::vector<int>& indices) const
{
    if (indices.empty())
        return topology.residues;

    std::vector<std::string> names;
    names.reserve(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        names.push_back(topology.residues[i]);
    }
    return names;
}

std::vector<int> Trajectory::get_residue_numbers (const std::vector<int>& indices) const
{
    if (indices.empty())
        return topology.residue_numbers;

    std::vector<int> numbers;
    numbers.reserve(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        numbers.push_back(topology.residue_numbers[i]);
    }
    return numbers;
}

std::vector<std::string> Trajectory::get_elements (const std::vector<int>& indices) const
{
    if (indices.empty())
        return topology.elements;

    std::vector<std::string> names;
    names.reserve(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        names.push_back(topology.elements[i]);
    }
    return names;
}

std::vector<int> Trajectory::get_chain_ids (const std::vector<int>& indices) const
{
    if (indices.empty())
        return topology.chain_ids;

    std::vector<int> chains;
    chains.reserve(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        chains.push_back(topology.chain_ids[i]);
    }
    return chains;
}

std::vector<std::vector<std::array<float,3>>> Trajectory::get_coordinates_stl (const std::vector<int>& indices) const
{
    std::vector<int> _indices;
    if (indices.empty()) {
        _indices = std::vector<int>(coordinates[0].size());
        int n = {0};
        std::generate(_indices.begin(), _indices.end(), [&n]{ return n++; });
    }

    std::vector<std::vector<std::array<float,3>>> coords;
    //NOT GOOD: TODO FIX!!! 
    for (int frame = 0; frame < n_frames; frame++) {
    //for (auto frame : _indices) {
        std::vector<std::array<float,3>> atom_coords;
        for (int atom = 0; atom < _indices.size(); atom++) {
            atom_coords.push_back({coordinates[frame][_indices[atom]][0], coordinates[frame][_indices[atom]][1], coordinates[frame][_indices[atom]][2]});
        }
        coords.push_back(atom_coords);
    }

    return coords;
}

//==============================================================================
//  Read different data formats
//==============================================================================

Topology Trajectory::read_PDB_topology (const std::string& topology_filename) const
{
    PDBReader pr;
    pr.read_file(topology_filename);

    Topology t;
    t.atom_names = pr.get_atom_names();
    t.elements = pr.get_elements();
    t.residues = pr.get_residue_names();
    t.residue_numbers = pr.get_residue_numbers();
    t.chain_ids = pr.get_chain_ids();

    return t;
}

Topology Trajectory::read_TOP_topology (const std::string& topology_filename) const
{
    TopReader tr;
    tr.read_file(topology_filename);

    Topology t;
    t.atom_names = tr.get_atom_names();
    t.elements = tr.get_elements();
    t.residues = tr.get_residue_names();
    t.residue_numbers = tr.get_residue_numbers();
    t.chain_ids = tr.get_chain_ids();

    return t;
}

std::vector<std::vector<std::array<float,3>>> Trajectory::read_netCDF_trajectory (const std::string& trajectory_filename) const
{
    TrajReader tr;
    tr.read_file(trajectory_filename);
    return tr.get_raw_trajectories();
}

std::vector<std::vector<std::array<float,3>>> Trajectory::read_PDB_trajectory (const std::string& trajectory_filename) const
{
    PDBReader pr;
    pr.read_file(trajectory_filename);
    std::vector<std::vector<std::array<float,3>>> coords;
    coords.push_back(pr.get_atom_coordinates());
    return coords;
}
// The PDBreader code is extracting all the atoms of the input PDB files 
// the code below is separating the different snapshots 
// this code is taking data as input which contains all the coordinates and splitting them and saving them
// in new vector array i.e. parsed atom. 
std::vector<std::vector<std::array<float,3>>> Trajectory::parse_raw_pdbtraj (int n_atoms, int n_dims, int n_frames, std::vector<std::vector<std::array<float,3>>> data)
{
    std::vector<std::vector<std::array<float,3>>> parsed_data;
    for (int frame = 0; frame < n_frames; frame++) {
        std::vector<std::array<float,3>> atom_coords;
        for (int atom = 0; atom < n_atoms; atom++) {
            std::array<float,3> coord = {
                data[0][n_atoms*frame + atom][0],
		data[0][n_atoms*frame + atom][1],
		data[0][n_atoms*frame + atom][2]
            };
            atom_coords.push_back(coord);
        }
        parsed_data.push_back(atom_coords);
    }

    return parsed_data;
} 
