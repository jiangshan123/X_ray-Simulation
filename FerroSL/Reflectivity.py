import imp  # python module to import other files
import matplotlib.pyplot as plt  # plotting
import math as m

# importing my own classes from a file
xpy = imp.load_source('xpy', r'/Users/Rui/Dropbox/Research/Stony Brook/Matt Dawber/Python Simulation/x-ray_simulator/xpy/src/xpy.py')

def degtoarc(deg):
    """converts degree to arcs
    deg: float
    returns: float
    """
    return deg / 180 * m.pi


# Def Data class
def arctodeg(arc):
    """converts arcs to degree
    arc: float
    returns: float
    """
    return arc * 180 / m.pi


def qz(theta2, Lambda=1.5406 * 10 ** -10):  # you can write powers in python as **| 3**2=9
    """calculates qz for a theta-2theta scan
    theta2: float (angle in degree)
    Lambda: float (wavelength in meters)
    """
    return 2 * m.pi / Lambda * m.sin(degtoarc(theta2) / 2) * 2  # q*sin(theta)*2


def l(q, a=3.905 * 10 ** -10):
    """calculates h from q and a
    a: float (substrate lattice parameter in meters)
    q: float (momentum transver in meters**-1)
    """
    return q * a / 2 / m.pi


def load_scans_tupple(path):
    """this function imports data from a simple format(2theta,intensety) from the file in 'path'
    path: string
    """
    data = []  # list to store our data
    data1 = []  # list to store our second data
    n = 0.0

    f = open(path, 'r')  # opens a file at the beginning
    line = f.readline()  # reads the next line of the file (until next \n)

    while (line) and n<5.0:  # this checks if there is something in line not NON
        words = line.split('  ')
        if words:  # just to be save that you could split the line
            data.append(float(words[0]))  # append adds one element to the end of the list
            data1.append(float(words[1]))  # append adds one element to the end of the list
            n = float(words[0])
        line = f.readline()

    return data, data1

def load_scans_tupple_2(path):
    """this function imports data from a simple format(2theta,intensety) from the file in 'path'
    path: string
    """
    data = []  # list to store our data
    data1 = []  # list to store our second data

    f = open(path, 'r')  # opens a file at the beginning
    line = f.readline()  # reads the next line of the file (until next \n)
    while (line):  # this checks if there is something in line not NON
        words = line.split(',')  # split the line into a list of strings for each ','
        if words:  # just to be save that you could split the line
            data.append(float(words[0]))  # append adds one element to the end of the list
            data1.append(float(words[1]))  # append adds one element to the end of the list
        line = f.readline()

    return data, data1


class X_Ray_Scan:  # everything that belongs to this class needs to be indented
    def __init__(self, path):
        self.import_data(path)

    def _calc_qz_(self):
        self.qz = [qz(theta2) for theta2 in self.twotheta]
        #self.qz = self.twotheta

    def _calc_l_(self):
        self.l = [l(qz) for qz in self.qz]

    def import_data(self, path):
        """this function imports data from a simple format(2theta,intensety) from the file in 'path'
        path: string
        """
        self.twotheta, self.intensety = load_scans_tupple(path)
        self._calc_qz_()
        self._calc_l_()

class X_Ray_Scan_2:  # everything that belongs to this class needs to be indented
    def __init__(self, path):
        self.import_data(path)

    def _calc_qz_(self):
        self.qz = [qz(theta2) for theta2 in self.twotheta]
        #self.qz = self.twotheta

    def _calc_l_(self):
        self.l = [l(qz) for qz in self.qz]

    def import_data(self, path):
        """this function imports data from a simple format(2theta,intensety) from the file in 'path'
        path: string
        """
        self.twotheta, self.intensety = load_scans_tupple_2(path)
        self._calc_qz_()
        self._calc_l_()

fig = plt.figure(figsize=(16,10), dpi=100)
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
axes.set_xlabel(r'$2\theta$' + ' [' + r'$^\circ$' + ']', fontsize=35)
axes.set_ylabel('int[A.U.]', fontsize=35)
axes.set_yscale("log")
axes.set_ylim(3*10**3,10**6)
axes.set_xlim(0.6,1.5)



# Experiment
path = r'/Users/Rui/Dropbox/Research/STONY BROOK/Matt Dawber/BTO_growthrate_Project/Raw Data/XRD/J0915161/J0915161_reflectivity_vs_10sec.UXD'
data = X_Ray_Scan_2(path)
axes.plot(data.twotheta, data.intensety, color="blue", label=r'Experimental Data',ls='-', marker='.', ms =17)

# Fitting
path2 = r'/Users/Rui/Dropbox/Research/STONY BROOK/Matt Dawber/BTO_growthrate_Project/Raw Data/XRD/J0915161/aa.UXD'
data2 = X_Ray_Scan(path2)
axes.plot(data2.twotheta, data2.intensety, color="red", label=r'Fitting',ls='-', marker='', ms =17,linewidth=3.0)

for tick in axes.xaxis.get_major_ticks():
    tick.label.set_fontsize(30)
for tick in axes.yaxis.get_major_ticks():
    tick.label.set_fontsize(30)
axes.grid(True)
axes.legend(loc = 1, fontsize=30, numpoints=1,markerscale=3)
#plt.title(r'$550^\circ C$ BTO on 3 Layer PTO Reflectivity Scan', fontsize=35, color="black")
#plt.show()
plt.savefig('/Users/Rui/Dropbox/Research/STONY BROOK/Matt Dawber/BTO_growthrate_Project/reflectivity_Simulation.pdf')