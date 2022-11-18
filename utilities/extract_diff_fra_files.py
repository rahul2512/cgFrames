import numpy as np
import sys

#Cytosine
#Thymine
#Thymine
#Adenine
#Guanine
#Thymine
#Thymine
#Cytosine
#Adenine
#Adenine
#Adenine
#Thymine
#Thymine
#Thymine
#Guanine
#Adenine
#Adenine
#Cytosine
#Thymine
#Adenine
#Adenine
#Guanine
#Cytosine
#Guanine
#
#Cytosine
#Thymine
#Thymine
#Adenine
#Guanine
#Thymine
#Thymine
#Cytosine
#Adenine
#Adenine
#Adenine
#Thymine
#Thymine
#Thymine
#Guanine
#Adenine
#Adenine
#Cytosine
#Thymine
#Adenine
#Adenine
#Guanine
#Cytosine

def main ():
    fra1 = np.genfromtxt(sys.argv[1], dtype=None, max_rows=48)
    fra2 = np.genfromtxt(sys.argv[2], dtype=None, max_rows=48)

    W_diff_x = [];
    W_diff_y = [];
    W_diff_z = [];
    C_diff_x = [];
    C_diff_y = [];
    C_diff_z = [];

    midframe_pos_fra1 = []
    for i in range(24):
        midframe_pos_fra1.append(np.array([0.5*(fra1[i][11]+fra1[i+24][11]), 0.5*(fra1[i][12]+fra1[i+24][12]), 0.5*(fra1[i][13]+fra1[i+24][13])]))
        #print midframe_pos[-1]

    midframe_pos_fra2 = []
    for i in range(24):
        midframe_pos_fra2.append(np.array([0.5*(fra2[i][11]+fra2[i+24][11]), 0.5*(fra2[i][12]+fra2[i+24][12]), 0.5*(fra2[i][13]+fra2[i+24][13])]))
        #print midframe_pos[-1]

    #for i in range(24):
    #    if fra1[i][0] == 1:
    #        W_diff_x.append((fra1[i][11]-midframe_pos_fra1[i][0])-(fra2[i][11]-midframe_pos_fra2[i][0]))
    #        W_diff_y.append((fra1[i][12]-midframe_pos_fra1[i][0])-(fra2[i][12]-midframe_pos_fra2[i][0]))
    #        W_diff_z.append((fra1[i][13]-midframe_pos_fra1[i][0])-(fra2[i][13]-midframe_pos_fra2[i][0]))
    #for i in range(24):
    #    if fra1[i][0] == 2:
    #        C_diff_x.append((fra1[i+24][11]-midframe_pos_fra1[i][0])-(fra2[i+24][11]-midframe_pos_fra2[i][0]))
    #        C_diff_y.append((fra1[i+24][12]-midframe_pos_fra1[i][0])-(fra2[i+24][12]-midframe_pos_fra2[i][0]))
    #        C_diff_z.append((fra1[i+24][13]-midframe_pos_fra1[i][0])-(fra2[i+24][13]-midframe_pos_fra2[i][0]))

    for i in range(24):
        if fra1[i][0] == 1:
            W_diff_x.append(fra1[i][11]-fra2[i][11])
            W_diff_y.append(fra1[i][12]-fra2[i][12])
            W_diff_z.append(fra1[i][13]-fra2[i][13])
    for i in range(24):
        if fra1[i][0] == 2:
            C_diff_x.append(fra1[i+24][11]-fra2[i+24][11])
            C_diff_y.append(fra1[i+24][12]-fra2[i+24][12])
            C_diff_z.append(fra1[i+24][13]-fra2[i+24][13])

    for i in range(24):
        print i, W_diff_x[i], W_diff_y[i], W_diff_z[i]
    print ' '
    for i in range(24):
        print i, C_diff_x[i], C_diff_y[i], C_diff_z[i]

if __name__ == "__main__":
    main()
