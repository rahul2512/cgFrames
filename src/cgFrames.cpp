#include "cgFrames.hpp"



template <class T>
bool is_set_equal (std::vector<T> set1, std::vector<T> set2)
{
    std::sort(set1.begin(), set1.end());
    std::sort(set2.begin(), set2.end());
    return std::includes(set1.begin(), set1.end(), set2.begin(), set2.end());
}

void cgFrames::read_ideal_bases (const std::string& filename)
{
    std::cout << "Reading ideal bases...";
    AtomGroupReader ag_reader;
    ag_reader.read_file(filename);
    ideal_bases = ag_reader.get_atom_groups();
    std::cout << "done" << std::endl; 
}

void cgFrames::read_trajectory (const std::string& trajectory_filename,
                                const std::string& topology_filename)
{
    trajectory.load(trajectory_filename, topology_filename);
}

void cgFrames::strip_trajectory ()
{
    std::vector<int> ib_atom_indices;
    for (atom_group ag : ideal_bases) {
        for (auto atom : ag.atoms) {
            auto atom_indices = trajectory.select_atom_with_name(atom.name);
            //append to other indices: (also copies duplicates atm
            ib_atom_indices.insert(ib_atom_indices.end(), atom_indices.begin(), atom_indices.end());
        }
    }
//    for (int j=0;j<ib_atom_indices.size();j++){
// 	std::cout <<  ib_atom_indices[j]<< '\n' << std::flush ; 
//    }	
    trajectory.remove(trajectory.select_inverse(ib_atom_indices));
}

std::vector<atom_group> cgFrames::find_groups_in_trajectory ()
{
    std::vector<atom_group> atom_groups;
    int group_indices = 0;

    std::cout << "Looking for atom groups in trajectory..." << std::flush;
    auto atom_names = trajectory.get_names();
    int amount_atoms = atom_names.size();
//    for (int i = 0; i < amount_atoms; i++) {
//	std::cout << i <<"     " << atom_names[i] <<"    "<< "\n" << std::flush ;
//	}
    for (int i = 0; i < amount_atoms; i++) {
        for (auto& ib : ideal_bases) {
            int len = ib.n_atoms;
	    //std::cout << i <<"     " <<len<<"     " << "\n" << std::flush ; 
            //if we go out of bounds, check the next ideal base
            if (i+len > amount_atoms)
		//std::cout << "halkat" << std::flush;
                continue;

            //construct set of atom names in ideal base:
            std::vector<std::string> ib_atom_names;
            ib_atom_names.reserve(15);
	    //std::cout << ib_atom_names[0] << "\n" << std::flush ;
            for (auto& atom : ib.atoms)
                ib_atom_names.push_back(atom.name);
//	        std::cout << i << "    " << len <<"   "<< atom_names[i] <<"\n" <<std::flush ;
            if (is_set_equal(std::vector<std::string>(atom_names.begin()+i, atom_names.begin()+i+len), ib_atom_names)) {
                atom_group group;
 		//std::cout << i <<"     " <<len<<"     " << atom_names[i] <<"    "<< amount_atoms<<"\n" << std::flush ;
                std::cout << "found something: ";
                std::cout << ib.name << std::endl;

                group.name = ib.name;
                group.index = group_indices;
                group_indices++;

                //get indices to put stuff in order:
                for (auto& ib_atom : ib.atoms) {
                    atom a;
                    for (int j = i; j < i+len; j++) {
                        if (atom_names.at(j) == ib_atom.name) {
                            a.name = atom_names.at(j);
                            a.md_index = j;
                        }
                    }
                    group.atoms.emplace_back(a);
                }
                atom_groups.emplace_back(group);
            }
        }
    }
    std::cout << "Done" << std::endl;
    std::cout << "Total number of groups found: " << atom_groups.size() << std::endl;

    return atom_groups;
}

std::vector<std::vector<atom_group>> cgFrames::get_atom_group_trajectories (const std::vector<int>& frame_indices)
{
    std::cout << "Reorganising data for fitting..." << std::flush;
    auto atom_groups_per_frame = find_groups_in_trajectory();
    std::vector<std::vector<atom_group>> atom_group_trajectories;

    int amount_frames = 0;
    if (frame_indices.size() != 0)
        amount_frames = frame_indices.size();
    else
        amount_frames = trajectory.get_number_of_md_frames();
    
    atom_group_trajectories.reserve(amount_frames);
    for (int i = 0; i < amount_frames; i++)
        atom_group_trajectories.emplace_back(atom_groups_per_frame);

    auto traj_coordinates = trajectory.get_coordinates_stl();
    int n_frames = traj_coordinates.size();

    //get atom coordinates:
    for (int frame_index = 0; frame_index < amount_frames; ++frame_index) {
        for (auto& group : atom_group_trajectories[frame_index]) {
            for (auto& atom : group.atoms) {
                atom.position = traj_coordinates[frame_index][atom.md_index];
            }
        }
    }
    std::cout << "Done" << std::endl;
    
    return atom_group_trajectories;
}

void cgFrames::fit_frames_to_trajectory (std::vector<std::vector<atom_group>>& ag_coords)
{
//#pragma omp parallel for
	//for (auto& md_frame : ag_coords) 
	for (auto ii = 0u; ii < ag_coords.size(); ++ii) 
	{
		auto& md_frame = ag_coords[ii];
		for (auto& group : md_frame) 
		{
			atom_group ideal_base;
			for (auto& ib : ideal_bases) 
			{
				if (ib.name == group.name) 
				{
					ideal_base = ib;
					break;
				}
			}
			fit.fit_base(group, ideal_base);
	//		std::cout << ii << "   ii   "<< std::flush;
		}
	}
}
