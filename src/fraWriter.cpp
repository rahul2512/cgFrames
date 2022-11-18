#include "fraWriter.hpp"

FraWriter::FraWriter () {}
FraWriter::~FraWriter () {}

//================================================================================
//  Setting outputs
//================================================================================

void FraWriter::set_output_names (const std::vector<std::pair<std::string,std::string>>& names)
{
    groups_file_names.clear();

    for (auto& name : names)
        groups_file_names.insert(name);
}

void FraWriter::add_output_name (const std::pair<std::string,std::string>& name)
{
    groups_file_names.insert(name);
}

void FraWriter::remove_output_name (const std::string& group_name)
{
    groups_file_names.erase(group_name);
}

void FraWriter::set_output_to_dna (const std::string& fra_file_name, const std::string& pfra_file_name)
{
    groups_file_names.clear();

    groups_file_names.insert({"Adenine",  fra_file_name});
    groups_file_names.insert({"Cytosine", fra_file_name});
    groups_file_names.insert({"Guanine",  fra_file_name});
    groups_file_names.insert({"Thymine",  fra_file_name});
    groups_file_names.insert({"Uracil",  fra_file_name});

    groups_file_names.insert({"Phosphate", pfra_file_name});
}

//================================================================================
//  Writing
//================================================================================

//When using the fraWriter, I KNOW that my data is DNA, and that bases in the same bp need to have the same index. Also know that strands begin and end with D*3/D*5,
//and that first strand is read in 5->3 direction, then the other one needs to be in 3->5 direction.
void FraWriter::write (std::vector<std::vector<atom_group>>& groups) const
{
    std::cout << "Writing frames...." << std::flush;
    for (int i = 0; i < groups.size(); ++i) {
        int number_of_groups = groups[i].size();
        //reverse the half of each set of groups:
        std::reverse(groups[i].begin()+number_of_groups/2, groups[i].end());
        int current_index = 1;
	int back_index = 2 ; 
        for (int j = 0; j < groups[i].size()/2 ; ++j) {
//	    std::cout << groups[i].size()/2 << std::flush;

            if (groups[i][j].name == "Phosphate")
                {groups[i][j].index = 1;
		groups[i][j].chain_index = back_index++;
		}
            else 
                {groups[i][j].index = current_index++; 
		groups[i][j].chain_index = 1;
                }
        }
        current_index = 1;
	back_index = 1;
        for (int j = groups[i].size()/2; j < groups[i].size(); ++j) {
            if (groups[i][j].name == "Phosphate")
                {groups[i][j].index = 2; 
		groups[i][j].chain_index = back_index++;
		}
            else 
                {groups[i][j].index = current_index++; 
		groups[i][j].chain_index = 2;            
		}
        }
    }


    //open file(s):
    std::unordered_map<std::string, std::ofstream> out_streams;
    for (int i = 0; i < groups[0].size(); ++i) {
        if (out_streams.find(groups[0][i].name) == out_streams.end()) {
            out_streams.insert(std::make_pair(groups[0][i].name, std::ofstream(groups_file_names.at(groups[0][i].name), std::ios::app)));
        }
    }

    for (int md_frame_index = 0; md_frame_index < groups.size(); ++md_frame_index) {

        for (auto& group : groups[md_frame_index]) {

            std::ofstream& file = out_streams.at(group.name);

            if (!out_streams.at(group.name).is_open()) {
                //first try to re-open.
                std::cout << "stream was closed" << std::endl;
                out_streams.at(group.name).open(groups_file_names.at(group.name), std::ios::app);
                //if still not open, nevermind
                if (!out_streams.at(group.name).is_open()) {
                    std::cout << "Opening output file failed. Exiting" << std::endl;
                    exit(1);
                }
            }
            //set stream:

            //std::cout << group.residue_number+1 << '\t' << group.chain_id+1 << '\t' <<group.residue_name << '\t' <<group.name << std::endl;
            //then write it:
            file << group.chain_index <<  "  " << group.index << "  ";
            //file << 1 << '\t' << group.index+1 << '\t';
            file << std::fixed << std::setprecision(7) << std::setw(13) << std::right;
            //file << std::scientific << std::setprecision(7) << std::setw(13) << std::right;
            file << group.R(0,0) << "  " << group.R(0,1) << "  " << group.R(0,2) << "  ";
            file << group.R(1,0) << "  " << group.R(1,1) << "  " << group.R(1,2) << "  ";
            file << group.R(2,0) << "  " << group.R(2,1) << "  " << group.R(2,2) << "  ";
            file << group.r(0)   << "  " << group.r(1)   << "  " << group.r(2) << '\n';
            file.flush();
        }
    }

    //close file(s)
    for (auto& s : out_streams)
        s.second.close();

    std::cout << "Done" << std::endl;
}
