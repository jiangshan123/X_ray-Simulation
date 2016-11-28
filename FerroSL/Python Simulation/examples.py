import imp #python module to import other files
%pylab inline #creates matplotlib plots inside the browser window
import matplotlib.pyplot as plt #plotting

#needed to use Latex inside the matplotlib figures
from matplotlib import rc
rc('font',**{'family':'sns-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#importing my own classes from a file
xpy = imp.load_source('xpy', 
		r'/home/lamandrionn/workspace/x-ray_simulator/xpy/src/xpy.py') 

# defining Experiment
l_min = 0.8 #first l point calculates
l_max = 1.05 #last l point calculates
l_step = 0.0002 #step sizes to go through the l range

wavelength = 1.5409E-10 #X-ray wavelength
direct = 2.0E-13 #can be used to adjust the intensety of the calculation to that of the data
background = 1 #number of background counts

########
# Film #
########

# defining the structures
STOc = xpy.CrystalStructure(path = xpy._SCRIPTPATH_+"Perovskites/STO_a.str") #loads a STO CrystalStructure 
PTO2c = xpy.create_structure(3.905, 4.156, 'Pb', 'Ti', 'Pb', 'Ti', 1, 0, 1)  #creates a PTO CrystalStructure with one unitcell
PTOthickc = xpy.create_structure(3.905, 4.156, 'Pb', 'Ti', 'Pb', 'Ti', 50, 0, 1) #creates a PTO CrystalStructure with 50 unitcell
STO = xpy.UnitCell(STOc) #puts the STO CrystalStructure into the Cell framework
PTO2 = xpy.UnitCell(PTO2c) #puts the PTO CrystalStructure into the Cell framework
PTOthick = xpy.UnitCell(PTOthickc) #puts the PTO CrystalStructure into the Cell framework

# defining the films

PTOThickstruck = xpy.Film(PTOthick,0,0,1) #creates the 50 unitcell film, note only one repetion!
PTOstruck = xpy.Film(PTO2,0,0,50) #creates the 50 layer film

#calculating the diffraction patterns
thickstruck = xpy.l_scan(xpy.ThickSample([substrate,PTOThickstruck]), STO, l_min = l_min, l_max = l_max, l_step = l_step, direct = direct, background = background, wavelength = wavelength)
struck =  xpy.l_scan(xpy.ThickSample([substrate,PTOstruck]), STO, l_min = l_min, l_max = l_max, l_step = l_step, direct = direct, background = background, wavelength = wavelength)

#make plot for the films
fig = plt.figure(figsize=(16,10), dpi=100) #creates the figure

#add a plot to the figure
axes = fig.add_axes([0., 0., 1, 1]) # left, bottom, width, height (range 0 to 1)

#add data sets to the figure
axes.plot([x[2] for x in thickstruck], [y[3] for y in thickstruck], color="blue", label='50 layer PTO unit cell',ls='*', marker='.', ms =15)
axes.plot([x[2] for x in struck], [y[3] for y in struck], color="green", label='50 layer PTO film',ls='*', marker='.', ms =15)

axes.set_xlabel(r'[$Q_z$]', fontsize=45) #define the x label
axes.set_ylabel('int[A.U.]', fontsize=45) #define the y label
axes.set_title('', fontsize=55) #define title
axes.set_yscale("log") #set the y scale to logarithmic plotting
axes.set_xlim(1.25*10**10,1.7*10**10) #set the x scale limits
axes.set_ylim(1,10**8) #set the y scale limits
axes.grid(True) #add grid lines to guide the eye
axes.legend(loc = 2, fontsize=40, numpoints=1,markerscale=3)# add the legend
#increase the tick label font size
for tick in axes.xaxis.get_major_ticks():
    tick.label.set_fontsize(40)
for tick in axes.yaxis.get_major_ticks():
    tick.label.set_fontsize(40)

#add text to the figure (here used to add the figure label)
axes.text(1.21*10**(10),8*10**7,'a', fontsize=65) #(x,y,text,fontsize) x,y position in plot coordinates!

#make plot for the film errors
fig = plt.figure(figsize=(16,10), dpi=100)
axes = fig.add_axes([0., 0., 1, 1]) # left, bottom, width, height (range 0 to 1)

axes.plot([x[2] for x in struck], [abs(float(a[3]) - float(b[3])) for a, b in zip(struck, thickstruck)], color="red", label='absolute error',ls='*', marker='.', ms =15)
axes.plot([x[2] for x in struck], [(abs(float(a[3]) - float(b[3])))/float(a[3]) for a, b in zip(struck, thickstruck)], label='relative error',ls='*', marker='.', ms =15)

axes.set_xlabel(r'[$Q_z$]', fontsize=45)
axes.set_ylabel('int[A.U.]', fontsize=45)
axes.set_yscale("log")
axes.text(1.202*10**(10),6*10**-8,'b', fontsize=65)
axes.legend(loc = 2, fontsize=40, numpoints=1,markerscale=3)
for tick in axes.xaxis.get_major_ticks():
    tick.label.set_fontsize(40)
for tick in axes.yaxis.get_major_ticks():
    tick.label.set_fontsize(40)
axes.grid(True)

################
# Superlattice #
################

# defining additional structures and films
PTOc = xpy.CrystalStructureCheck(path = xpy._SCRIPTPATH_+"Perovskites/PTO.str")
PTO = xpy.UnitCell(PTOc)
STOfilm = xpy.Film(STO,0,0,6)
PTOfilm = xpy.Film(PTO,0,0,2)
Bilayer = xpy.ThickSample([PTOfilm,STOfilm])
superfilm = xpy.Film(Bilayer, 0, 0 ,15)
supersample = xpy.ThickSample([PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm,PTOfilm,STOfilm])
Superfilm = xpy.ThickSample([substrate, superfilm])
Supersample = xpy.ThickSample([substrate, supersample])

#calculating the diffraction patterns
SuperFilm = xpy.l_scan(Superfilm, STO, l_min = l_min, l_max = l_max, l_step = l_step, direct = direct, background = background, wavelength = wavelength)
SuperSample = xpy.l_scan(Supersample, STO, l_min = l_min, l_max = l_max, l_step = l_step, direct = direct, background = background, wavelength = wavelength)

#make plot for the superlattice
fig = plt.figure(figsize=(16,10), dpi=100)
axes = fig.add_axes([0., 0., 1, 1]) # left, bottom, width, height (range 0 to 1)
axes.plot([x[2] for x in SuperFilm], [y[3] for y in SuperFilm], color="blue", label='superlattice created using the Film class',ls='*', marker='.', ms =15)
axes.plot([x[2] for x in SuperSample], [y[3] for y in SuperSample], color="green", label='superlattice created using the ThickSample class',ls='*', marker='.', ms =15)
axes.set_xlabel(r'[$Q_z$]', fontsize=45)
axes.set_ylabel('int[A.U.]', fontsize=45)
axes.set_yscale("log")
axes.set_ylim(1,10**8)
axes.text(1.21*10**(10),8*10**7,'a', fontsize=65)
axes.legend(loc = 2, fontsize=40, numpoints=1,markerscale=3)
for tick in axes.xaxis.get_major_ticks():
    tick.label.set_fontsize(40)
for tick in axes.yaxis.get_major_ticks():
    tick.label.set_fontsize(40)
axes.grid(True)

#make plot for the superlattice errors
fig = plt.figure(figsize=(16,10), dpi=100)
axes = fig.add_axes([0., 0., 1, 1]) # left, bottom, width, height (range 0 to 1)
axes.plot([x[2] for x in SuperFilm], [abs(float(a[3]) - float(b[3])) for a, b in zip(SuperFilm, SuperSample)], color="red", label='absolute error',ls='*', marker='.', ms =15)
axes.plot([x[2] for x in SuperFilm], [(abs(float(a[3]) - float(b[3])))/float(a[3]) for a, b in zip(SuperFilm, SuperSample)], label='relative error',ls='*', marker='.', ms =15)
axes.set_xlabel(r'[$Q_z$]', fontsize=45)
axes.set_ylabel('int[A.U.]', fontsize=45)
axes.set_yscale("log")
axes.text(1.202*10**(10),6*10**-8,'b', fontsize=65)
axes.legend(loc = 2, fontsize=40, numpoints=1,markerscale=3)
for tick in axes.xaxis.get_major_ticks():
    tick.label.set_fontsize(40)
for tick in axes.yaxis.get_major_ticks():
    tick.label.set_fontsize(40)
axes.grid(True)

