#include "TrajReader.hpp"

template <class T>
bool is_set_equal (std::vector<T> set1, std::vector<T> set2)
{
    std::sort(set1.begin(), set1.end());
    std::sort(set2.begin(), set2.end());
    return std::includes(set1.begin(), set1.end(), set2.begin(), set2.end());
}

void TrajReader::read_file (const std::string& filename)
{
    std::cout << "Reading data from netCDF file..." << std::flush;
    try {
        netCDF::NcFile data_file = netCDF::NcFile(filename, netCDF::NcFile::read);
        netCDF::NcVar v = data_file.getVar("coordinates");

        if(v.isNull()) {
            std::cout << "not reading correctly"<< std::endl;
            exit(1);
        }
        //Get some metadata to determine the number of frames, and atoms and such.
        int n_params = v.getDimCount();
        int n_dims = 0;
        int n_atoms = 0;
        int n_frames = 0;
        for (int d = 0; d < n_params; d++) {
            if (v.getDim(d).getName() == "frame")
                n_frames = v.getDim(d).getSize();
            if (v.getDim(d).getName() == "atom")
                n_atoms = v.getDim(d).getSize();
            if (v.getDim(d).getName() == "spatial")
                n_dims = v.getDim(d).getSize();
        }

        amount_dims = n_dims;
        amount_atoms = n_atoms;
        amount_frames = n_frames;

        float* data = new float[n_frames*n_atoms*n_dims];

        //Get the data:
        v.getVar(data);
        raw_trajectories = parse_raw_traj(n_atoms, n_dims, n_frames, data);

        //std::cout << "raw traj info: " << std::endl;
        //std::cout << "raw traj frames: " << std::endl;
        //std::cout << raw_trajectories.size() << std::endl;
        //std::cout << "raw traj atoms: " << std::endl;
        //std::cout << raw_trajectories.at(0).size() << std::endl;
        //std::cout << raw_trajectories[0][0][0] << '\t' << raw_trajectories[0][0][1] << '\t' << raw_trajectories[0][0][2] << std::endl;

        delete [] data;
    } catch (netCDF::exceptions::NcException& e) {
        std::cout << e.what() << std::endl;
        exit(1);
    }
    std::cout << "Done" << std::endl;
}

std::vector<std::vector<std::array<float,3>>> TrajReader::parse_raw_traj (int n_atoms, int n_dims, int n_frames, float* data)
{
    std::vector<std::vector<std::array<float,3>>> parsed_data;
    for (int frame = 0; frame < n_frames; frame++) {
        std::vector<std::array<float,3>> atom_coords;
        for (int atom = 0; atom < n_atoms; atom++) {
            std::array<float,3> coord = {
                data[n_dims*n_atoms*frame + n_dims*atom + 0],
                data[n_dims*n_atoms*frame + n_dims*atom + 1],
                data[n_dims*n_atoms*frame + n_dims*atom + 2]
            };
            atom_coords.push_back(coord);
        }
        parsed_data.push_back(atom_coords);
    }

    return parsed_data;
}

void TrajReader::dispose_of_hydrogen (const std::vector<std::string>& atom_names)
{
    //first find every hydrogen atom:
    std::vector<int> indices;
    for (int i = 0; i < atom_names.size(); i++) {
        if (atom_names[i][0] == 'H')
            indices.push_back(i);
    }
    //then remove them from the data:
    for (int f = 0; f < amount_frames; f++) {
        for (int i = indices.size()-1; i >= 0; i--) {
            raw_trajectories[f].erase(raw_trajectories[f].begin()+indices.at(i));
        }
    }
}

void TrajReader::dispose_of_ions (const std::vector<std::string>& atom_names)
{
    //first find every hydrogen atom:
    std::vector<int> indices;
    for (int i = 0; i < atom_names.size(); i++) {
        if (atom_names[i] == "K+" || atom_names[i] == "Na+" || atom_names[i] == "Cl-" )
            indices.push_back(i);
    }
    //then remove them from the data:
    for (int f = 0; f < amount_frames; f++) {
        for (int i = indices.size()-1; i >= 0; i--) {
            raw_trajectories[f].erase(raw_trajectories[f].begin()+indices.at(i));
        }
    }
}

std::vector<std::vector<std::array<float,3>>> TrajReader::get_raw_trajectories ()
{
    return raw_trajectories;
}
