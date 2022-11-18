import numpy as np
import sys

def main ():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    seqfile = ""
    position_seqfile = 0

    fra = np.genfromtxt(infile, dtype=None, max_rows=48)
    #seq = ['G' for i in range(max(zip(*fra)[1]))]
    seq = ['G','C','T','T','A','G','T','T','C','A','A','A','T','T','T','G','A','A','C','T','A','A','G','C']
    print seq
    #fra = np.genfromtxt(infile, delimiter=5, dtype=int)
    bp = parse_fra(zip(*fra), seq)
    json = make_json(bp)
    write_json(outfile, json);


def parse_fra(fra, seq):
    chain = [i for i in fra[0]]
    index = [i for i in fra[1]]
    R00 = [float(i)/1.0000 for i in fra[2]]
    R01 = [float(i)/1.0000 for i in fra[3]]
    R02 = [float(i)/1.0000 for i in fra[4]]
    R10 = [float(i)/1.0000 for i in fra[5]]
    R11 = [float(i)/1.0000 for i in fra[6]]
    R12 = [float(i)/1.0000 for i in fra[7]]
    R20 = [float(i)/1.0000 for i in fra[8]]
    R21 = [float(i)/1.0000 for i in fra[9]]
    R22 = [float(i)/1.0000 for i in fra[10]]
    x = [float(i)/1.0000 for i in fra[11]]
    y = [float(i)/1.0000 for i in fra[12]]
    z = [float(i)/1.0000 for i in fra[13]]

    basepairs = []
    nbp = len(index)/2
    for i in range(nbp):
        bp = {'R' : [R00[i], R01[i], R02[i], R10[i], R11[i], R12[i], R20[i], R21[i], R22[i]],
                'r' : [x[i],y[i],z[i]],
                'Rc' : [R00[i+nbp], R01[i+nbp], R02[i+nbp], -R10[i+nbp], -R11[i+nbp], -R12[i+nbp], -R20[i+nbp], -R21[i+nbp], -R22[i+nbp]],
                'rc' : [x[i+nbp], y[i+nbp], z[i+nbp]], 'sequence' : seq[i]}
        basepairs.append(bp)
    return basepairs

def make_json (bp):
    json = '{\n\
            \"Data\" : [ {\n\
            \t\t\"basepair\" : ['
    for i in range(len(bp)):
        json += '\t{ \n\
                \t\t\t\t\"R\" : ' + str(bp[i]['R']) + ',\n\
                \t\t\t\t\"r\" : ' + str(bp[i]['r']) + ',\n\
                \t\t\t\t\"Rc\" : ' + str(bp[i]['Rc']) + ',\n\
                \t\t\t\t\"rc\" : ' + str(bp[i]['rc']) + ',\n\
                \t\t\t\t\"sequence\" : ' + '\"' + bp[i]['sequence'] + '\"' + '\n\
                \t\t\t}'
        if i < (len(bp)-1):
            json += ', '
    json += '\t\t\t], \"energy\" : 0.0\n\
            \t\t}\n\
            \t]\n\
            }'

    return json

def write_json (outfile,json):
    with open(outfile, 'w') as out:
        out.write(json)

if __name__ == '__main__':
    main()

