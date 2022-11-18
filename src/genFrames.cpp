#include <iostream>
#include "cgFrames.hpp"

int main (int argc, char** argv)
{
    //takes as input:
    //  -T  trajectory file
    //  -t  topology file
    //  -i  ideal bases file
    //  -o  output file(s)

    //FOR NOW:
    //  trajectory file
    //  topology file
    //  ideal bases file
    //  fra output filename
    //  pfra output filename
    //get all of these:
    if (argc < 6) {
        std::cout << "Usage: " << argv[0] << " <trajectory filename> <topology filename> <ideal bases filename> <fra out filename> <pfra out filename> <fitting algo (optional)>" << std::endl;
        exit(1);
    }

    std::string traj_file = argv[1];
    std::string top_file = argv[2];
    std::string ib_file = argv[3];
    std::string fra_out_file = argv[4];
    std::string pfra_out_file = argv[5];
//    std::string fit_algo = argv[6];

    cgFrames cgf;

    cgf.read_trajectory(traj_file, top_file);
    cgf.read_ideal_bases(ib_file);

    //gets rid of all atoms that are not in the ideal bases file
    cgf.strip_trajectory();


    auto gat = cgf.get_atom_group_trajectories();

    //get frames
    cgf.fit_frames_to_trajectory(gat);

    cgf.output.set_output_to_dna(fra_out_file, pfra_out_file);
    cgf.output.write(gat);

    return 0;
}
