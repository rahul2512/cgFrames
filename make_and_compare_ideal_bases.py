import numpy as np
import subprocess
import sys
import glob
import time

def generate_ideal_bases_file (folder, Ax, Ay, Cx, Cy, Gx, Gy, Tx, Ty):
    filename = folder+"ideal_bases-A_"+str(Ax)+"_"+str(Ay)+"-C"+str(Cx)+"_"+str(Cy)+"-G"+str(Gx)+"_"+str(Gy)+"-T"+str(Tx)+"_"+str(Ty)+".txt"
    with open(filename, "w") as f:
        f.write("#entries are as follows:\n"
                "#<group name> <# atoms in group> [ideal frame: 9x orientation + 3x position]\n"
                "# <x> <y> <z> <atom name>\n"
                "#\n"
                "Adenine 11 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 "+str(Ax)+" "+str(Ay)+" 0.0000\n"
                "-2.497963270852318 5.389504179143331 -5.551115123125783e-17 C1'\n"
                "-1.286577654496827 4.522557013344271 -5.551115123125783e-17 N9\n"
                "0.041448512680411 4.886554049508940 0.001998788600148 C8\n"
                "0.862626669580093 3.883338036610646 0.001331109569129 N7\n"
                "0.033335728421553 2.771051874622240 4.886954884535666e-04 C5\n"
                "0.298311935152744 1.392142742959572 2.291118707407475e-04 C6\n"
                "1.610999331226802 0.929985147423594 0 N6\n"
                "-0.758766225646914 0.558232846858952 -0.001302764712450 N1\n"
                "-1.985886075239647 1.077036918897750 -0.001839114848758 C2\n"
                "-2.353923620065637 2.338772962882393 -0.001796589074074 N3\n"
                "-1.277547971696661 3.152750831970960 -4.163336342344337e-17 C4\n"
                "Thymine 10 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 "+str(Tx)+" "+str(Ty)+" 0.0000\n"
                "-2.49797385251954 5.38951175207988 -8.32667268468867e-17 C1'\n"
                "-1.28657765449683 4.52255701334427 -9.71445146547012e-17 N1\n"
                "-0.0308318929617097 5.06966398538589 0.00106227587351299 C6\n"
                "1.06705884597297 4.29544362667478 0.00182968625221110 C5\n"
                "2.46409441406427 4.97805677636900 0.00100000000000000 C7\n"
                "0.954415605430206 2.85551051188829 0.00131997916448753 C4\n"
                "1.93792690492813 2.13682539678559 0 O4\n"
                "-1.49052504046889 3.16358559151261 -1.38777878078145e-16 C2\n"
                "-2.56735108282079 2.63243335990271 0 O2\n"
                "-0.345449648812713 2.39168814612310 0.00125651287566437 N3\n"
                "Guanine 12 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 "+str(Gx)+" "+str(Gy)+" 0.0000\n"
                "-1.28657765449683 4.52255701334427 -6.93889390390723e-17 N9\n"
                "0.0426423566174657 4.88330289243608 0.00130968377641127 C8\n"
                "0.871086535669193 3.86785583503164 0.00131167498235545 N7\n"
                "0.0311247947646614 2.75498381388369 0.000526268903276811 C5\n"
                "0.347544708728062 1.37173063868634 0.000999083717970736 C6\n"
                "1.55199075361480 0.923057827003115 0 O6\n"
                "-0.811580042335409 0.582460960590514 -0.000109613967279027 N1\n"
                "-2.10468805144495 1.06601640189687 4.02518859410478e-05 C2\n"
                "-2.95201095238939 0.112605094074547 -0.00100000000000000 N2\n"
                "-2.39995031196363 2.36417953688744 -0.000724381917298830 N3\n"
                "-1.28738844874338 3.14401190725862 -6.93889390390723e-17 C4\n"
                "-2.49796842535900 5.38950786804709 -1.11022302462516e-16 C1'\n"
                "Cytosine 9 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 "+str(Cx)+" "+str(Cy)+" 0.0000\n"
                "-1.28657765449683 4.52255701334427 2.77555756156289e-17 N1\n"
                "-1.47665054038848 3.14316712685502 0 C2\n"
                "-2.62254356529968 2.68442396996229 0.00100000000000000 O2\n"
                "-0.385057702897589 2.33540381348435 0.000345992810207529 N3\n"
                "0.848496113805832 2.85477978961241 0.00158616112557461 C4\n"
                "1.88303685950519 2.01968489662960 0.00100000000000000 N4\n"
                "1.06531043633039 4.27112225105204 0.00155664258340361 C5\n"
                "-0.0382837938843514 5.06162143324159 0.00103501592415939 C6\n"
                "-2.49796842535900 5.38950786804710 -5.55111512312578e-17 C1'\n"
                "Phosphate 5 0.28880532 -0.50008344 0.81639941 -0.40811277 0.70707284 0.57748763 -0.8659639 -0.50010651 0.00 0.0 0.0 0.0\n"
                "0.000 0.000 0.000 P\n"
                "1.518 0.000 -0.537 O3'\n"
                "-0.759 -1.315 -0.537 O5'\n"
                "-0.698 1.208 -0.493 OP1\n"
                "0.000 0.000 1.480 OP2")
        return filename

def gen_all_ib_files():
    for Ax in np.arange(-3.0,3.0,1.0):
        for Ay in np.arange(-3.0,3.0,1.0):
            for Cx in np.arange(-3.0,3.0,1.0):
                for Cy in np.arange(-3.0,3.0,1.0):
                    for Gx in np.arange(-3.0,3.0,1.0):
                        for Gy in np.arange(-3.0,3.0,1.0):
                            for Tx in np.arange(-3.0,3.0,1.0):
                                for Ty in np.arange(-3.0,3.0,1.0):
                                    name = generate_ideal_bases_file("./test_files/ideal_bases/", Ax, Ay, Cx, Cy, Gx, Gy, Tx, Ty)
                                    #subprocess.Popen(["./genFrames", "./test_files/Palin_141_1.1.ions.nc", "./test_files/Palin_141_1.ions.top", name, "./test_files/ideal_bases_out_fra/out-A_"+str(Ax)+"_"+str(Ay)+"-C"+str(Cx)+"_"+str(Cy)+"-G"+str(Gx)+"_"+str(Gy)+"-T"+str(Tx)+"_"+str(Ty)+".fra"])

def run_all():
    files = glob.glob("test_files/ideal_bases/ideal_bases*.txt")
    outputs = ["./test_files/ideal_bases_out_fra/"+f[0:-4].split('/')[-1] for f in files]
    #print files
    #print outputs
    for i in range(0,len(files)/50, 11):
        children = []
        for j in range(0,11):
            children.append(subprocess.Popen(("./genFrames", "./test_files/Palin_141_1.1.ions.nc", "./test_files/Palin_141_1.ions.top", files[i], outputs[i]+".fra", outputs[i]+".pfra")))
        for c in children:
            c.wait()



if __name__ == "__main__":
    gen_all_ib_files()
    #run_all()

