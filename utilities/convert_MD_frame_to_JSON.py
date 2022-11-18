import numpy as np
import sys
import mdtraj as md

"""
Usage: python convert_MD_frame_to_JSON.py <trajectory file> <topology file> <frame no.> <out_file>
"""

def main ():
    traj = md.load(sys.argv[1], top=sys.argv[2])
    frame_no = int(sys.argv[3])
    out_file = sys.argv[4]

    #print traj
    #print traj.xyz
    #print traj[frame_no].xyz[0]
    #print traj.top.atom(1).name
    #print traj.top.atom(1).residue.name
    #print traj.top.residue(1).name
    #coords = traj[frame_no].xyz
    traj.remove_solvent()
    for i in range(traj.n_atoms):
        print "\t"+str(i)+"\t"+traj.top.atom(i).name+"\t"+traj.top.atom(i).residue.name
    save_as_JSON(traj[frame_no].xyz[0], traj.topology, out_file)


def save_as_JSON (MD_frame, topology, outfile):
    out = ""
    out += "{\n"
    out += "\t\"atoms\" : [\n"
    skipped = 0
    for i in range(len(MD_frame)):
        if (topology.atom(i).name == 'K+'
            or topology.atom(i).name == 'Cl-'
            or topology.atom(i).name == "C2'"
            or topology.atom(i).name == "C3'"
            or topology.atom(i).name == "C4'"
            or topology.atom(i).name == "O4'"
            or topology.atom(i).name == "C5'"
            or topology.atom(i).name[0] == 'H'):
            skipped += 1
            continue
        else:
            out += "\t\t{\n"
            out += "\t\t\t\"name\" : \""+topology.atom(i).name+"\",\n"
            out += "\t\t\t\"index\" : "+str(i-skipped)+",\n"
            out += "\t\t\t\"residue\" : \""+topology.atom(i).residue.name+"\",\n"
            out += "\t\t\t\"coords\" : ["+str(MD_frame[i][0])+", "+str(MD_frame[i][1])+", "+str(MD_frame[i][2])+"]\n"
        if i != (len(MD_frame)-skipped-1):
            out += "\t\t},\n"
        else:
            out += "\t\t}\n"
    out += "\t]\n"
    out += "}\n"

    with open(outfile, 'w') as f:
        f.write(out)



if __name__ == "__main__":
    main()
