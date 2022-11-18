import numpy as np

A = {'C1' : np.array([  1.57340, -2.41044, -0.12190]), 'N9'  : np.array([ 0.37786, -1.52945, -0.00531]), 'C4' : np.array([0.38110, -0.16006, -0.04011])}

T = {'C1' : np.array([1.95363, -1.45659, -0.19073]), 'N1' : np.array([0.75787, -0.57587, -0.07419]), 'C2' : np.array([0.97215, 0.78019, -0.13405])}

G = {'C1' : np.array([1.58195, -2.39594, -0.12320]), 'N9' : np.array([0.37714, -1.52785, -0.00520]), 'C4' : np.array([0.37551, -0.14975, -0.04020])}

C = {'C1' : np.array([1.94866, -1.45161, -0.19029]), 'N1' : np.array([0.74385, -0.58352, -0.07229]), 'C2' : np.array([0.93020, 0.79520, -0.12929 ])}



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
        vec = nuc['N9'] + d*np.dot(rot_t1, (nuc['N9']-nuc['C1']))#/np.linalg.norm(nuc['N9']-nuc['C1']))
    elif type == "pyrimidine":
        vec = nuc['N1'] + d*np.dot(rot_t1, (nuc['N1']-nuc['C1']))#/np.linalg.norm(nuc['N1']-nuc['C1']))

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

    print 'A:', A_bR
    print 'T:', T_bR
    print 'G:', G_bR
    print 'C:', C_bR
