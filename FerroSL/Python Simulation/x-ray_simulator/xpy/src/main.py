import math as m
import xpy
import xpyc
import cProfile

l_min = 0.8
l_max = 1.05
l_step = 0.000005

h_min = 1-0.07
h_max = 1+0.07
h_step = 0.0007

wavelength = 1.23984*10**-10
direct = 2.0E-12
background = 0

STOc = xpy.CrystalStructureCheck(path = xpy._SCRIPTPATH_+"Perovskites/STO.str")
PTOc = xpy.CrystalStructureCheck(path = xpy._SCRIPTPATH_+"Perovskites/PTOs.str")
STO = xpy.UnitCell(STOc)
PTO = xpy.UnitCell(PTOc)
substrate = xpy.Film(STO)
film = xpy.Film(PTO,0,0,50)
#film = xpy.SuperCell(PTOc,1,1,50)

experiment = xpy.Experiment(substrate, STOc, direct = direct, background = background, wavelength = wavelength)

sample = xpy.ThickSample([film])

experiment = xpy.Experiment(sample, STO, direct = direct, background = background, wavelength = wavelength)

for i in experiment.l_scan(l_min, l_max, l_step):
    print(i)