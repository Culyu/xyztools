# -*- coding: utf-8 -*-
"""
Bondorientation calculator, histogram version

To run:
python bondorientation.py xyzfile.xyz vectorfile kind cutoff binsize

direction = x, y, z, xy

X = all kinds

@author: sampokulju
"""
# Importing
import sys
import tools
import numpy as np
import math

# Get the filenames from arguments
xyzfile = sys.argv[1]
vectorfile = sys.argv[2]

# Get the atom kinds and cutoff
center = sys.argv[3]
cutoff = float(sys.argv[4])


# Make data object
xyz = tools.dataread(xyzfile)
vectors = tools.dataread(vectorfile)

# Make data array like "a list in a list"
xyz.xyzfile()
vectors.parsecol()

# Make bins
binsize = float(sys.argv[5])
binnum = int(math.ceil(float(vectors.data[2][2])/binsize))
xybins = []
zbins = []
for num in range(binnum):
    binvalue = round(num * binsize, 2)
    xybins.append([binvalue, 0])
    zbins.append([binvalue, 0])

"""
Make neighborlist
"""

print("Making neighborlist...")

mult = [0, 1]

neighlist = []
atomindex = 0
ccount = 0
zsum = 0
for atom in xyz.data:
    kind = atom[0]
    x = float(atom[1])
    y = float(atom[2])
    z = float(atom[3])
    atomlist = [[kind, x, y, z]]
    if(kind == 'C'):
        zsum = zsum+z
        ccount = ccount+1
    neighbors = []
    for multx in mult:
        for multy in mult:
            for multz in mult:
                # print "Multiplier: ",multx,', ', multy,', ',multz
                neighindex = 0
                for neigh in xyz.data:
                    kindneigh = neigh[0]
                    xneigh = (float(neigh[1]) +
                              multx*float(vectors.data[0][0]) +
                              multy*float(vectors.data[1][0]) +
                              multz*float(vectors.data[2][0]))
                    yneigh = (float(neigh[2]) +
                              multx*float(vectors.data[0][1]) +
                              multy*float(vectors.data[1][1]) +
                              multz*float(vectors.data[2][1]))
                    zneigh = (float(neigh[3]) +
                              multx*float(vectors.data[0][2]) +
                              multy*float(vectors.data[1][2]) +
                              multz*float(vectors.data[2][2]))
                    distvec = [x-xneigh, y-yneigh, z-zneigh]
                    distance = np.sqrt(distvec[0]**2 +
                                       distvec[1]**2 +
                                       distvec[2]**2)
                    if(distance < cutoff):
                        if(atomindex != neighindex and kindneigh != 'C'):
                            edge = abs(multx) + abs(multy) + abs(multz)
                            neighinfo = [kindneigh, xneigh, yneigh, zneigh, edge]
                            neighbors.append(neighinfo)
                    neighindex += 1
    atomlist.append(neighbors)
    neighlist.append(atomlist)
    atomindex += 1

if(ccount != 0):
    graphh = zsum/ccount

print("Calculating bond orientations and placing atoms to bins...")

"""
Go trough neighborlist and calculate orientation
"""

for nn in neighlist:
    centerkind = nn[0][0]
    if (centerkind == center or center == 'X' and centerkind != 'C'):
        centerx = nn[0][1]
        centery = nn[0][2]
        centerz = nn[0][3]
        nnlist = nn[1]
        for bond in nnlist:
            bondvector = [bond[1]-centerx, bond[2]-centery, bond[3]-centerz]
            bondcenter = [centerx + bondvector[0]/2,
                          centery + bondvector[1]/2,
                          centerz + bondvector[2]/2]
            proz = abs(bondvector[2])
            proxy = np.sqrt(bondvector[0]**2 + bondvector[1]**2)
            bondlength = np.sqrt(bondvector[0]**2 +
                                 bondvector[1]**2 +
                                 bondvector[2]**2)
            if(bond[4] != 0):
                edgescale = 2.0
            else:
                edgescale = 1.0

            if(ccount != 0):
                shiftedz = bondcenter[2]-graphh
            else:
                shiftedz = bondcenter[2]
            # Calculate bin
            placeto = int(math.floor(shiftedz/binsize))
            if(proz >= proxy):
                binsum = zbins[placeto][1]
                newbinsum = float(binsum) + 0.5*edgescale
                zbins[placeto][1] = newbinsum
            else:
                binsum = xybins[placeto][1]
                newbinsum = float(binsum) + 0.5*edgescale
                xybins[placeto][1] = newbinsum

# Print plotting data
filename = "xyorientation.dat"
f = open(filename, 'w')
try:
    f.write("#Bond orientation respect to XY \n")
    for k in range(len(xybins)):
        for value in xybins[k]:
            value = str(value)
            f.write("%s   " % value)
        f.write("\n")
finally:
    f.close()

filename = "zorientation.dat"
f = open(filename, 'w')
try:
    f.write("#Bond orientation respect to Z \n")
    for k in range(len(zbins)):
        for value in zbins[k]:
            value = str(value)
            f.write("%s   " % value)
        f.write("\n")
finally:
    f.close()
