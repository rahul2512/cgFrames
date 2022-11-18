import numpy as np


A = {'C1' : np.array([-0.2479, 0.5346, 0.000]),'N9' : np.array([-0.1291, 0.4498, 0.000]),'C4' : np.array([-0.1267, 0.3124, 0.000])}

T = {'C1' : np.array([-0.2481, 0.5354, 0.000]), 'N1' : np.array([-0.1284, 0.45, 0.000]), 'C2' : np.array([-0.1462, 0.3135, 0.000])}

G = {'C1' : np.array([-0.2477, 0.5399, 0.000]), 'N9' : np.array([-0.1289, 0.4551, 0.000]), 'C4' : np.array([-0.1265, 0.3177, 0.000])}

C = {'C1' : np.array([-0.2477, 0.5402, 0.000]), 'N1' : np.array([-0.1285, 0.4542, 0.000]), 'C2' : np.array([-0.1472, 0.3158, 0.000])}


"""
A = {'C1' : np.array([-2.497963270852318, 5.389504179143331, -5.551115123125783e-17]),'N9' : np.array([-1.286577654496827, 4.522557013344271, -5.551115123125783e-17]),'C4' : np.array([-1.277547971696661, 3.152750831970960, -4.163336342344337e-17])}

T = {'C1' : np.array([-2.49797385251954, 5.38951175207988, -8.32667268468867e-17]), 'N1' : np.array([-1.28657765449683, 4.52255701334427, -9.71445146547012e-17]), 'C2' : np.array([-1.49052504046889, 3.16358559151261, -1.38777878078145e-16])}

G = {'C1' : np.array([-2.49796842535900, 5.38950786804709, -1.11022302462516e-16]), 'N9' : np.array([-1.28657765449683, 4.52255701334427, -6.93889390390723e-17]), 'C4' : np.array([-1.28738844874338, 3.14401190725862, -6.93889390390723e-17])}

C = {'C1' : np.array([-2.49796842535900, 5.38950786804710, -5.55111512312578e-17]), 'N1' : np.array([-1.28657765449683, 4.52255701334427, 2.77555756156289e-17]), 'C2' : np.array([-1.47665054038848, 3.14316712685502, 0.00000])}
"""

def bN (nuc, type):
    normal = 0;
    if type == "purine":
        normal = np.cross((nuc['N9']-nuc['C1']),(nuc['N9']-nuc['C4']))
    elif type == "pyrimidine":
        normal = np.cross((nuc['N1']-nuc['C1']),(nuc['N1']-nuc['C2']))

    norm2 = np.linalg.norm(normal)

    return normal/norm2

def bL (nuc, type):
    t2 = (-54.41/180.0)*np.pi
    rot_t2 = np.array([[np.cos(t2), -np.sin(t2),0], [np.sin(t2), np.cos(t2),0],[0,0,1]])
    vec = 0
    if type == "purine":
        vec = np.dot(rot_t2, (nuc['N9']-nuc['C1']))
    elif type == "pyrimidine":
        vec = np.dot(rot_t2, (nuc['N1']-nuc['C1']))
    norm2 = np.linalg.norm(vec)

    return vec/norm2

def bR (nuc, type):
    t1 = (141.47/180.0)*np.pi
    d = 4.702 #angstrom
    rot_t1 = np.array([[np.cos(t1), -np.sin(t1),0], [np.sin(t1), np.cos(t1),0],[0,0,1]])
    vec = 0
    if type == "purine":
        vec = d*np.dot(rot_t1, (nuc['N9']-nuc['C1'])/np.linalg.norm(nuc['N9']-nuc['C1']))
    elif type == "pyrimidine":
        vec = d*np.dot(rot_t1, (nuc['N1']-nuc['C1'])/np.linalg.norm(nuc['N1']-nuc['C1']))

    return vec

if __name__ == "__main__":
    A_bN = bN(A, 'purine')
    A_bL = bL(A, 'purine')
    A_bD = np.cross(A_bN, A_bL)
    A_bR = bR(A, 'purine')

    G_bN = bN(G, 'purine')
    G_bL = bL(G, 'purine')
    G_bD = np.cross(G_bN, G_bL)
    G_bR = bR(G, 'purine')

    T_bN = bN(T, 'pyrimidine')
    T_bL = bL(T, 'pyrimidine')
    T_bD = np.cross(T_bN, T_bL)
    T_bR = bR(T, 'pyrimidine')

    C_bN = bN(C, 'pyrimidine')
    C_bL = bL(C, 'pyrimidine')
    C_bD = np.cross(C_bN, C_bL)
    C_bR = bR(C, 'pyrimidine')

    #careful about order of vectors?
    A_R = np.array([A_bD, A_bL, A_bN])
    G_R = np.array([G_bD, G_bL, G_bN])
    T_R = np.array([T_bD, T_bL, T_bN])
    C_R = np.array([C_bD, C_bL, C_bN])

    print 'A:'
    print '\tbR:',A_bR
    print 'R:',A_R
    print 'T:'
    print 'R:',T_R
    print 'G:'
    print 'R:',G_R
    print 'C:'
    print 'R:',C_R


    #print np.dot(A_R,np.transpose(A_R))
    #print np.dot(G_R,np.transpose(G_R))
    #print np.dot(T_R,np.transpose(T_R))
    #print np.dot(C_R,np.transpose(C_R))
