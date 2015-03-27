# -*- coding: utf-8 -*-
"""
Script to add partial (Bader) charges to XYZ-file for visualization.

To run:
python bader.py xyzfile.xyz baderfile.dat

baderfile is ACF.dat
"""

# Importing
import sys
import tools


# Get the filenames from arguments
xyzfile = sys.argv[1]
baderfile = sys.argv[2]

# Make data objects
xyz = tools.dataread(xyzfile)
bader = tools.dataread(baderfile)

# Split lines to columns
xyz.parsecol()
bader.parsecol()

# Check the number of atoms in data
atomnum = int(xyz.data[0][0])

# Delete headers from both data
xyzhead = 2
xyz.delhead(xyzhead)

baderhead = 2
bader.delhead(baderhead)

# Delete footers or extra lines from data
baderfoot = 4
bader.delfoot(baderfoot)

xyzextra = len(xyz.data)-atomnum
xyz.delfoot(xyzextra)

# Add column from Bader-data to xyz

xyzbader = xyz
xyzbader.addcolumn(bader.data, 4)

# Add number of atoms and comment line to make correct XYZ-filetype

xyzbader.xyzheader()

# Print the file

filename = "xyzbader.xyz"
xyzbader.dataprint(filename)
