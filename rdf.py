# -*- coding: utf-8 -*-
"""
Radial distribution function (RDF) calculator

To run:
python rdf.py xyzfile.xyz vectorfile binsize

vectorfile:
CELLVEC1X CELLVEC1Y CELLVEC1Z
CELLVEC2X CELLVEC2Y CELLVEC2Z
CELLVEC3X CELLVEC3Y CELLVEC3Z

@author: sampokulju
"""
# Importing
import sys
import tools
import numpy as np
import itertools
import math

# Get the filenames from arguments
xyzfile = sys.argv[1]
vectorfile = sys.argv[2]

# Get the bin size
binsize = float(sys.argv[3])


# Make data object
xyz = tools.dataread(xyzfile)
vectors = tools.dataread(vectorfile)

# Make data array like "a list in a list"
xyz.xyzfile()
vectors.parsecol()

# Cutoff = the half of the shortest cell vector
vec1length = np.sqrt(float(vectors.data[0][0])**2 +
                     float(vectors.data[0][1])**2 +
                     float(vectors.data[0][2])**2)
vec2length = np.sqrt(float(vectors.data[1][0])**2 +
                     float(vectors.data[1][1])**2 +
                     float(vectors.data[1][2])**2)
vec3length = np.sqrt(float(vectors.data[2][0])**2 +
                     float(vectors.data[2][1])**2 +
                     float(vectors.data[2][2])**2)

minlength = vec1length
if(vec2length < minlength):
    minlength = vec2length

if(vec3length < minlength):
    minlength = vec3length

maxrad = minlength/2

# Number of bins
binfloat = maxrad/binsize
binnum = math.trunc(binfloat)
cutoff = binnum*binsize

rbins = np.linspace(binsize, cutoff, num=binnum)


"""
Find all kinds, count how many and make pairs (Ge-Ge, Ge-Te...so on)
"""
kinds = []
for atomline in xyz.data:
    kind = atomline[0]
    if(kind not in kinds and kind != 'C'):
        kinds.append(kind)

pairs = list(itertools.product(kinds, kinds))

kindnumbers = []
for kind in kinds:
    kindsum = sum(1 for atomline in xyz.data if atomline[0] == kind)
    kindnumbers.append([kind, kindsum])

print(kindnumbers)

"""
Make RDF weight
"""
mult = [-1, 0, 1]

partialcoord = []

for pair in pairs:
    binsdata = []
    for value in rbins:
        binsdata.append([value, 0])

    atomindex = 0
    center = pair[0]
    compare = pair[1]
    print('Sort neighbors to bins for ' + center + ' - ' + compare)
    # Check number of atoms in bins
    for atom in xyz.data:
        kind = atom[0]
        x = float(atom[1])
        y = float(atom[2])
        z = float(atom[3])
        if(kind == center):
            # Check periodic images
            for multx in mult:
                for multy in mult:
                    for multz in mult:
                        neighindex = 0
                        for neigh in xyz.data:
                            kindneigh = neigh[0]
                            xneigh = (float(neigh[1]) +
                                      multx * float(vectors.data[0][0]) +
                                      multy * float(vectors.data[1][0]) +
                                      multz * float(vectors.data[2][0]))
                            yneigh = (float(neigh[2]) +
                                      multx * float(vectors.data[0][1]) +
                                      multy * float(vectors.data[1][1]) +
                                      multz * float(vectors.data[2][1]))
                            zneigh = (float(neigh[3]) +
                                      multx * float(vectors.data[0][2]) +
                                      multy * float(vectors.data[1][2]) +
                                      multz * float(vectors.data[2][2]))
                            distvec = [x-xneigh, y-yneigh, z-zneigh]
                            distance = np.sqrt(distvec[0]**2 + distvec[1]**2 +
                                               distvec[2]**2)
                            if(distance < cutoff):
                                if(atomindex != neighindex and kindneigh != 'C'
                                   and kindneigh == compare):
                                    binnotfound = True
                                    index = 0
                                    while(binnotfound):
                                        binmax = binsdata[index][0]
                                        if(distance <= binmax):
                                            numberinbin = binsdata[index][1]+1
                                            binsdata[index][1] = numberinbin
                                            binnotfound = False
                                        index += 1
                            neighindex += 1
        atomindex += 1

    for kind in kindnumbers:
        element = kind[0]
        if(element == compare):
            comparenum = kind[1]
        if(element == center):
            centernum = kind[1]

    # Calculate RDF
    print('Calculate RDF')
    rdf = []
    for atomsnum in binsdata:
        ninbin = atomsnum[1]
        rmax = float(atomsnum[0])
        rmin = rmax - binsize
        average = float(ninbin)/float(centernum)
        weight = (average/(4.0/3.0 * math.pi * ((float(rmax))**3 -
                  (float(rmin))**3))) / (float(comparenum) / (vec1length * vec2length * vec3length))
        rdf.append([rmax, weight])

    # Calculate integral to first min
    integ = 0
    for line in rdf:
        radius = float(line[0])
        weight = line[1]
        if(radius <= 3.21):
            integ = (integ + weight * 4.0/3.0 * math.pi * ((radius)**3 -
                     (radius-binsize)**3) * (float(comparenum) / (vec1length * vec2length * vec3length)))

    # Just two decimals
    integ = round(integ, 2)
    partialcoord.append([center, compare, integ])

    print('Writing data file')

    # Print RDF plotting data

    filename = center + "-" + compare + "-rdf.dat"
    f = open(filename, 'w')
    f.write("#RDF weight for " + center + " - " + compare + "  \n")
    for k in range(len(rdf)):
        for value in rdf[k]:
            value = str(value)
            f.write("%s   " % value)
        f.write("\n")


"""
Print partial coord data
"""

geoname = xyzfile[:-4]
filename = geoname + "-partialcoord.dat"
f = open(filename, 'w')
try:
    f.write("#Partial coordination numbers for  " + xyzfile + "  \n")
    for k in range(len(partialcoord)):
        for value in partialcoord[k]:
            value = str(value)
            f.write("%s   " % value)
        f.write("\n")
finally:
    f.close()
