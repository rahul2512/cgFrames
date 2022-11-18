import numpy as np
import sys

def main ():
    fra1 = np.genfromtxt(sys.argv[1], dtype=None, max_rows=48)
    fra2 = np.genfromtxt(sys.argv[2], dtype=None, max_rows=48)

    diff = 0;
    for i in range(48):
        diff += (fra1[i][11]-fra2[i][11])**2 + (fra1[i][12]-fra2[i][12])**2 + (fra1[i][13]-fra2[i][13])**2

    print diff

if __name__ == "__main__":
    main()
