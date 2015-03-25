class dataread:

    """
    Class for reading data from ASCII files and parsing and storing
    it to a list in a list to form array like structure with lines and
    columns.

    Also writing the same data in xyz format back to file.
    """

    # Read the data from files to list
    def __init__(self, filename):
        self.data = []
        self.data = open(filename, 'r').readlines()

    # Delete the header lines, number of lines is defined in main program
    def delhead(self, head):
        head = int(head)
        del self.data[0:head]

    # Delete the footer line, number of lines is defined in main program
    def delfoot(self, foot):
        if(foot != 0):
            foot = int(foot)
            del self.data[-foot:]

    # Split the lines containing one string to lists to form "columns"
    def parsecol(self):
        self.data = [line.split() for line in self.data]

    # Add one column from another data list to list in this object
    def addcolumn(self, datalist, column):
        for i in range(len(datalist)):
            line = datalist[i]
            addcell = line[column]
            self.data[i].append(addcell)

    # Add XYZ filetype header with number of atoms and empty commment line
    def xyzheader(self):
        atomsnum = len(self.data)
        self.data.insert(0, [" "])
        self.data.insert(0, [atomsnum])

    # Print the datalist to file
    def dataprint(self, filename):
        f = open(filename, 'w')
        for i in range(len(self.data)):
            for word in self.data[i]:
                word = str(word)
                f.write("%s   " % word)
            f.write("\n")

    # Add one line to file
    def addline(self, line):
        self.data.append(line)

    # Parse the XYZ
    def xyzfile(self):
        # Split lines to columns
        self.parsecol()

        # Check the number of atoms in data
        atomsnum = int(self.data[0][0])

        # Delete header from data
        xyzhead = 2
        self.delhead(xyzhead)

        # Delete possible extra lines from data
        extra = len(self.data)-atomsnum
        self.delfoot(extra)

    # Column modifier
    def coladdmod(self, datalist, column, height):
        for i in range(len(datalist)):
            line = datalist[i]
            addcell = float(line[column])
            addcell = (addcell-height)*(-5)
            self.data[i].append(addcell)

    # Sort according to column
    def numfloat(self, colnum):
        for i in range(len(self.data)):
            line = self.data[i]
            self.data[i] = line[:1] + [float(j) for j in line[1:]]


class datacollector:

    """
    To collect data from Bader data and calculate the sum of one kind
    """
    # Create the object with the list for kinds
    def __init__(self):
        self.kinds = []

    # Make a list for new kind and and data there
    def newkind(self, element, data):
        self.kinds.append(element)
        self.__dict__[element] = []
        data = float(data)
        self.__dict__[element].append(data)

    # Add data to existing list if kind is already found or make new
    # with newkind method if new kind is found
    def adddata(self, element, data):
        if element in self.kinds:
            data = float(data)
            self.__dict__[element].append(data)
        else:
            self.newkind(element, data)


class cubedata:

    """
    To collect data from cube-filetype to plot with i.e. gnuplot
    """
    def __init__(self, filename, zmin, zmax):
        self.data = []
        self.numberofatoms = 0
        self.cellvectors = []
        self.voxeln = []
        self.parsedata(filename, zmin, zmax)

    # New line from x y z values
    def addline(self, x, y, z, voxvalue):
        newline = [x, y, z, voxvalue]
        self.data.append(newline)

    # Print the datalist to file
    def dataprint(self, filename):
        f = open(filename, 'w')
        for i in range(len(self.data)):
            for word in self.data[i]:
                word = str(word)
                f.write("%s   " % word)
            f.write("\n")

    # Parse out the number of atoms and cell size
    def parsedata(self, filename, zmin, zmax):
        with open(filename) as f:

            # Line number
            ln = 0

            # Check metadata
            for commentline in f:
                # Number of atoms
                if ln == 2:
                    words = [commentline.split()]
                    self.numberofatoms = int(words[0][0])
                    print self.numberofatoms
                # Cell vectors
                if ln == 3 or ln == 4 or ln == 5:
                    words = [commentline.split()]
                    vnum = int(words[0][0])
                    self.voxeln.append(vnum)
                    vector = [float(words[0][1])*0.529177,
                              float(words[0][2])*0.529177,
                              float(words[0][3])*0.529177]
                    self.cellvectors.append(vector)
                if ln == 5:
                    break

                ln += 1
            # Max voxel numbers
            numx = self.voxeln[0]
            numy = self.voxeln[1]
            numz = self.voxeln[2]

            xcount = 0
            ycount = 0
            zcount = 0

            ln = 1

            for line in f:
                # Parse voxel data
                if ln > self.numberofatoms:
                    numbers = [line.split()]

                    for j in range(len(numbers[0])):
                        voxvalue = float(numbers[0][j])

                        dx = (xcount * self.cellvectors[0][0] +
                              ycount * self.cellvectors[1][0] +
                              zcount * self.cellvectors[2][0])
                        dy = (xcount * self.cellvectors[0][1] +
                              ycount * self.cellvectors[1][1] +
                              zcount * self.cellvectors[2][1])
                        dz = (xcount * self.cellvectors[0][2] +
                              ycount * self.cellvectors[1][2] +
                              zcount * self.cellvectors[2][2])

                        if dz >= zmin and dz <= zmax:
                                self.addline(dx, dy, dz, voxvalue)
                        zcount += 1

                        # Print zcount
                        if zcount == numz:
                            zcount = 0
                            ycount += 1

                        if ycount == numy:
                            ycount = 0
                            xcount += 1

                ln += 1

        # Symmetrical, copy data to form original
        xmult = numx / xcount - 1
        if xmult > 1:
            for row in range(len(self.data)):
                xcoord = self.data[row][0]
                ycoord = self.data[row][1]
                zcoord = self.data[row][2]
                vox = self.data[row][3]
                i = 1
                while i <= xmult:
                    xadd = xcoord + xcount * self.cellvectors[0][0] * i
                    self.addline(xadd, ycoord, zcoord, vox)
                    i += 1
