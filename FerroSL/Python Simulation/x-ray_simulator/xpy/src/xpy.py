import math as m
import cmath as mc
import fractions
import random
import multiprocessing as multi 

words = __file__.split("/")
_SCRIPTPATH_ = ""
for i in range(len(words)):
    if ( i == len(words)-3): break
    _SCRIPTPATH_ += words[i]+"/"
    
_NUMBER_ = (int, float) # np is not working in 3.3 yet ,np.floating, np.integer)



########
# Atom #
########

class Atom(object):
    
    def __init__(self, form, pos = None , N = 1, B = 0):
        """initializes the class data and checks data types
          
            pos: 3 element list of relative atom position default: [0.0 ...]  
            struc: Structure Factor of the Atom default: N/A
            B: is the Debye-Waller factor (Temperature dependent movement)
            N: probability that the atom is present in the a unit cell
        """
        if (pos):
            self.pos = pos
        else:
            self.pos = [0.0,0.0,0.0] 
        self.form = form
        self.B = B
        self.N = N
        #check data types
        self._Atom__type_check__()    
        
    def __str__(self):
        """Returns the data in a structured form.
           
           returns:
           atom: atom x=pos[1] y=pos[2] z= pos[3] 
        """
        
        stri = ("atom: " + self.form.atom + " x=" + str(self.pos[0]) +
                " y=" + str(self.pos[1]) + " z=" + str(self.pos[2]) + 
                " B=" + str(self.B) + " N=" + str(self.N)
                )
        return stri
        
    def __type_check__(self):
        """checks the class data to have the right types
        
            pos    -> list
            pos[i] -> Number
            struc  -> StructureFactor
            B -> Number
            N -> Number
        """
        assert isinstance(self.form, FormFactor) , 'form needs to be a FormFactor'
        assert (isinstance(self.pos, list) or (len(self.pos) != 3)) , 'a needs to be a list with 3 elements'
        for i in range(len(self.pos)):
            assert (isinstance(self.pos[i], _NUMBER_)) , 'pos' + str(i) + 'needs to be a number'
        assert (isinstance(self.B , _NUMBER_)) , 'B' + str(i) + 'needs to be a number'
        assert (isinstance(self.N , _NUMBER_)) , 'N' + str(i) + 'needs to be a number'
            
    _Atom__type_check__ = __type_check__ #private copy of the function to avoid overload

##############
# FormFactor #
##############

class FormFactor(object): 
    
    def __init__(self,atom="",a=None,b=None, i=1.0,q=0.0,value=0.0,path=None):
        """initializes the class data and checks data types
        
            atom:  chemical Symbol                   default: ""
            q:     absolute value of the wave vector default: 0.0   
            value: last calculated value             default: 0.0   
            a: 5 element list of prefactors          default: [0.0 ...]  
            b: 4 element list of exponents           default: [0.0 ...]
        """
        self.q = q
        self.value = value
        self.atom = atom
        self.i = i
        if (a):
            self.a = a
        else:
            self.a = [0.0,0.0,0.0,0.0,0.0,0.0]
        if(b): 
            self.b = b
        else:
            self.b = [0.0,0.0,0.0,0.0,0.0]
        #load data from path overwrites given data        
        if path:
            self.load_file(path)
        #check data types
        self._FormFactor__type_check__()
        
    def __str__(self):
        """Returns the data in a structured form.
        
            returns:
            element: 'atom'
            a1 = 'a[0]; a2 = 'a[1]; a3 = 'a[2]; a4 = 'a[3]; a5 = 'a[4] a6 = a[5];
            b1 = 'b[0]; b2 = 'b[1]; b3 = 'b[2]; b4 = 'b[3] b5 = b[4];
        """
        stri = "element:" + self.atom + "\n"
        for i in range(len(self.a)): 
            stri += 'a' + str(i+1) + " = " + str(self.a[i]) + "; "
        stri += "\n"
        for i in range(len(self.b)): 
            stri += 'b' + str(i+1) + " = " + str(self.b[i]) + "; "
        stri += "\n i: " + str(self.i)
        return stri
        
    def __type_check__(self):
        """checks the class data to have the right types.
        
            q -> float
            value -> float
            atom -> str
            a -> list
            a[i] -> float
            b -> list
            b[i] -> float
            i -> float
        """
        assert isinstance(self.q, _NUMBER_) , 'q needs to be a Number'
        assert isinstance(self.value, _NUMBER_) , 'value needs to be a Number'
        assert isinstance(self.atom, str) , 'atom needs to be a string'
        assert (isinstance(self.a, list) or (len(self.a) != 6)) , 'a needs to be a list with 5 elements'
        for i in range(len(self.a)):
            assert (isinstance(self.a[i], _NUMBER_)) , 'a' + str(i) + 'needs to be a Number'
        assert (isinstance(self.b, list) or (len(self.b) != 5)) , 'b needs to be a list with 4 elements' 
        for i in range(len(self.b)):
            assert (isinstance(self.b[i], _NUMBER_)) , 'b' + str(i) + 'needs to be a Number'
        assert isinstance(self.i, _NUMBER_) , 'q needs to be a Number'

    _FormFactor__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def get_value(self, q):
        """returns the structure factor value for the given q
        
           input q: absolute value of the wave vector^2 /16/m.pi/m.pi*10**-20     
        """
        #assert isinstance(q,self.q.__class__), 'q needs to be ' + str(self.q.__class__)
        if (self.q != q): 
            self.value =(self.a[0] * m.exp(-self.b[0] * q) + 
                         self.a[1] * m.exp(-self.b[1] * q) + 
                         self.a[2] * m.exp(-self.b[2] * q) + 
                         self.a[3] * m.exp(-self.b[3] * q) +
                         self.a[4] * m.exp(-self.b[4] * q) + 
                         self.a[5] + 1j * self.i
                        )
            self.q = q
        return self.value
        
    def load_file(self, path):
        """loads the structure factor the file in path"""
        f = open(path, 'r')
        line = f.readline()
        while (line):
            words = line.split()
            if words:
                if words[0] == "atom": self.atom = words[1]
                elif words[0] == "a1": self.a[0]= float(words[1])
                elif words[0] == "a2": self.a[1]= float(words[1])
                elif words[0] == "a3": self.a[2]= float(words[1])
                elif words[0] == "a4": self.a[3]= float(words[1])
                elif words[0] == "a5": self.a[4]= float(words[1])
                elif words[0] == "c": self.a[5]= float(words[1])
                elif words[0] == "b1": self.b[0]= float(words[1])
                elif words[0] == "b2": self.b[1]= float(words[1])
                elif words[0] == "b3": self.b[2]= float(words[1])
                elif words[0] == "b4": self.b[3]= float(words[1])
                elif words[0] == "b5": self.b[4]= float(words[1])
                elif words[0] == "i": self.i = float(words[1])
            line = f.readline()
        f.close()

####################
# CrystalStructure #
####################

class CrystalStructure(object):
    
    def __init__(self, lattice_parameters = None, atoms = None, damping=5.98E+4, path = None ):
        """initializes the class data and checks data types
          
            lattice_parameters: 3 element list default: [0.0 ...]  
            atoms: list of Atom instances default:[]
            damping: damping coefficient of the structure
            q: scattering vector default: [0.0,0.0,0.0]
            path: path to a file with structure information default: None
        """
        self.damping = damping
        self._factor_ = 10**-20/16/m.pi/m.pi
        if (atoms):
            self.atoms = atoms
        else:
            self.atoms = []
        if (lattice_parameters):
            self._lattice_parameters_ = lattice_parameters
        else:
            self._lattice_parameters_ = [1.0,1.0,1.0]
        self._structure_factor_ = 0.0+0.0j
        self._q_ = [0.0,0.0,0.0]
        #load data from path overwrites given data        
        if path:
            self.load_file(path)
        #check data types
        self._CrystalStructure__type_check__()
             
    def __str__(self):
        """Returns the data in a structured form.
        
        """
        stri = ("Lattice parameters: a=" + str(self._lattice_parameters_[0]) +
                " b="+ str(self._lattice_parameters_[1]) +
                " c="+ str(self._lattice_parameters_[2]) + "\n"
                )
        stri += "damping: " + str(self.damping) + '\n'
        for atom in self.atoms:
            stri += atom.__str__() + "\n"
        return stri
    
    def __type_check__(self):
        """checks the class data to have the right types
                  
            lattice_parameters -> list
            lattice_parameters[i] -> float
            atoms -> list
            atoms[i] -> Atom
            damping -> float
            _wavelength_ -> float
            _q_ -> list
            _q_[i[ -> float
            _formfactor_ -> complex
        """
        assert isinstance(self.damping, _NUMBER_) , 'damping needs to be a Number'
        assert (isinstance(self._lattice_parameters_, list) or (len(self._lattice_parameters_) != 3)) , 'lattice_parameters needs to be a list with 3 elements' 
        for i in range(len(self._lattice_parameters_)):
            assert (isinstance(self._lattice_parameters_[i], _NUMBER_)) , 'lattice_parameters[' + str(i) + '] needs to be a number'
        assert (isinstance(self.atoms, list)) , 'atoms needs to be a list ' 
        for i in range(len(self.atoms)):
            assert (isinstance(self.atoms[i], Atom)) , 'atoms[' + str(i) + '] needs to be a Atom'
    
    _CrystalStructure__type_check__ = __type_check__ #private copy of the function to avoid overload
            
    def _add_atom_(self, usedfile, atomtype):
        """adds the next atom in the file to the atom list """
        form = FormFactor(path=_SCRIPTPATH_+"/atoms/" + atomtype + ".at")
        pos = [0.0,0.0,0.0];
        line = usedfile.readline()
        count = 0
        B = 0
        N = 1
        while (line):
            words = line.split()
            if words:
                if words[0] == "B":
                    global B 
                    B = float(words[1])
                elif words[0] == "N":
                    global N 
                    N = float(words[1])
                elif words[0]   == "x": 
                    pos[0] = float(words[1])
                    count += 1
                elif words[0] == "y": 
                    pos[1] = float(words[1])
                    count += 1
                elif words[0] == "z": 
                    pos[2] = float(words[1])
                    count += 1
                if (count == 3): break 
            line = usedfile.readline()
        newatom = Atom(form, pos, N, B)
        self.add_atom(newatom)    

    def add_atom(self, atom):
        """adds atom to the atoms list"""
        assert isinstance(atom, Atom) , 'atom needs to be an Atom'
        for oldatom in self.atoms:
            if(atom.form.atom == oldatom.form.atom): 
                atom.form = oldatom.form
                break
        self.atoms.append(atom)

    def load_file(self, path):
        """loads the structure from the file in path"""
        f = open(path, 'r')
        line = f.readline()
        while (line):
            words = line.split()
            if words:
                if words[0]   == "element": self._add_atom_(f, words[1])
                elif words[0] == "a": self._lattice_parameters_[0] = float(words[1])
                elif words[0] == "b": self._lattice_parameters_[1] = float(words[1])
                elif words[0] == "c": self._lattice_parameters_[2] = float(words[1])
                elif words[0] == "damping": self.damping = float(words[1])
                
            line = f.readline()
        f.close()
    
    def save_file(self, path):
        """saves the structure to the file path"""
        f = open(path, 'w')
        f.write('a ' + str(self._lattice_parameters_[0]) + '\n')
        f.write('b ' + str(self._lattice_parameters_[1]) + '\n')
        f.write('c ' + str(self._lattice_parameters_[2]) + '\n')
        for atom in self.atoms:
            f.write('element ' + str(atom.form.atom) + '\n')
            f.write('B ' + str(atom.B) + '\n')
            f.write('N ' + str(atom.N) + '\n')
            f.write('x ' + str(atom.pos[0]) + '\n')
            f.write('y ' + str(atom.pos[1]) + '\n')
            f.write('z ' + str(atom.pos[2]) + '\n')
        f.close()
        
    def _calc_structure_factor_(self, q, sinomega):
        """calculates the structure factor for the structure 
        
            q: scattering vector
            sinomega: sin(angle between incident beam and sample surface
        """
        if (self.atoms == []):
            return 0.0
        structure_factor = 0.0 + 0.0j
        mq = sum([x**2 for x in q])
        q2 = mq * self._factor_
        damp = self.damping / sinomega * self._lattice_parameters_[2] #damping
        for atom in self.atoms:
            multi = 0.0
            for i in range(3):
                multi += -1j * (1 - atom.pos[i]) * self._lattice_parameters_[i] * q[i] 
            structure_factor +=(atom.form.get_value(q2) * atom.N *
                                mc.exp(-atom.B * q2 + multi - damp * (1 - atom.pos[2]))     
                               )
        return structure_factor
                  
    def get_structure_factor(self, q, sinomega):
        """calculates the structure factor for the structure 
        
            q: scattering vector
            sinomega: sin(angle between incident beam and sample surface
        """
        return self._calc_structure_factor_(q, sinomega)
    
    @property
    def a(self):
        """a lattice parameter"""
        return self._lattice_parameters_[0]

    @a.setter
    def a(self, value):
        self._lattice_parameters_[0] = value
        
    @property
    def b(self):
        """b lattice parameter"""
        return self._lattice_parameters_[1]

    @b.setter
    def b(self, value):
        self._lattice_parameters_[1] = value
        
    @property
    def c(self):
        """c lattice parameter"""
        return self._lattice_parameters_[2]

    @c.setter
    def c(self, value):
        self._lattice_parameters_[2] = value

#########################
# CrystalStructureCheck #
#########################

class CrystalStructureCheck(CrystalStructure):
    
    def get_structure_factor(self, q, sinomega):
        """calculates the structure factor for the structure 
        
            q: scattering vector
            sinomega: sin(angle between incident beam and sample surface
        """
        #assert isinstance(q,self._q_.__class__), 'q needs to be ' + str(self._q_.__class__)
        if (self._q_ != q or self._sinomega_ != sinomega): 
            self._structure_factor_ = self._calc_structure_factor_(q, sinomega)
            self._q_ = q[:]
            self._sinomega_ = sinomega
        return self._structure_factor_

########
# Cell #
########

class Cell(object):
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the UnitCell.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return self._calc_reflection_(q, sinomega)
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the UnitCell.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        return self._calc_phase_shift_(q, sinomega)
    
    def get_width(self):
        """returns the width of the UnitCell."""
        return self._calc_width_()
    
    def get_depth(self):
        """returns the depth of the UnitCell"""
        return self._calc_depth_()
    
    def get_thickness(self):
        """returns the thickness of the UnitCell."""
        return self._calc_thickness_() 
    
    @property
    def a(self):
        """a lattice parameter"""
        return self.get_width()
        
    @property
    def b(self):
        """b lattice parameter"""
        return self.get_depth()
        
    @property
    def c(self):
        """c lattice parameter"""
        return self.get_thickness()
    
    @property
    def damping(self):
        """c lattice parameter"""
        return self._calc_damping_()


############
# UnitCell #
############

class UnitCell(Cell):
    
    def __init__(self, crystal):
        """initializes the class data and checks data types
          
           crystal: CrystalStructure of the film  
        """
        self.crystal= crystal
        self._icq_ = self._ibq_ = self._iaq_ = 0.0
        self._eiaq_ = self._eibq_ = self._eicq_ = 0.0
        self._damping_ = 0.0
        self._amp_ =  self._calc_amplitude_()
    
    def __type_check__(self):
        """checks the class data to have the right types
        
           crystal -> Crystalstructure 
        """
        assert isinstance(self.crystal, CrystalStructure) , 'crystal needs to be a Crystalstructure'
        assert isinstance(self.eps, _NUMBER_) , 'eps needs to be a number'
        
    _UnitCell__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def __str__(self):
        """Returns the data in a structured form.
           
           returns its crystal structure 
        """
        return self.crystal.__str__()
    
    def _calc_amplitude_(self):
        """returns calculates the scattering amplitude"""
        R = 2.1879e-15 #R= 2.1879e-15 is the classical electron radius
        amp1d = 4 * m.pi * R / self.a / self.b 
        return amp1d 
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the UnitCell."""        
        r = self._amp_ * self.crystal.get_structure_factor(q, sinomega)
        return r
    
    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the UnitCell
        
            q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        self._set_icq_(q, sinomega)
        return [self._eiaq_, 
                self._eibq_, 
                self._eicq_ * self._damping_
               ]
    
    def _calc_damping_(self):
        """returns the damping factor of one unitCell of the underlying cristal"""
        return self.crystal.damping
     
    def _calc_width_(self):
        """returns the width of the UnitCell"""
        return self.crystal.a
    
    def _calc_depth_(self):
        """returns the breath of the UnitCell"""
        return self.crystal.b
    
    def _calc_thickness_(self):
        """returns the thickness of the UnitCell"""
        return self.crystal.c
    
    def _set_icq_(self,  q, sinomega):
        """recalculates the values for i*(latticeparameters*q*n)"""
        factor = self.damping * self.crystal.c / sinomega
        self._iaq_ = -1J * self.crystal.a * q[0]
        self._ibq_ = -1J * self.crystal.b * q[1] 
        self._icq_ = -1J * self.crystal.c * q[2]
        self._eiaq_ = mc.exp(self._iaq_)
        self._eibq_ = mc.exp(self._ibq_)
        self._eicq_ = mc.exp(self._icq_)
        self._damping_ = m.exp(-factor)
        
#############
# SuperCell #
#############

class SuperCell(UnitCell):

    def __init__(self, crystal, wide=0, deep=0, high=0):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           layers:  number of layers in the film
        """
        
        self._high_ = high
        self._deep_ = deep
        self._wide_ = wide
        super(self.__class__, self).__init__(crystal)
        self.eps = 0.01   
        self._SuperCell__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           wide -> int
           deep -> int
           high -> int 
        """
        assert isinstance(self._high_, (int)) , 'high needs to be a int'
        assert isinstance(self._deep_, (int)) , 'deep needs to be a int'
        assert isinstance(self._wide_, (int)) , 'wide needs to be a int'
        
    _SuperCell__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _calc_direction_(self, eiaq, wide):
        if wide == 0:
            if (eiaq.real > 1 - self.eps):
                return 1.0
            return 0.0
        return self._calc_out_direction(eiaq, wide, 0.9999999999)
    
    def _calc_out_direction(self, eiaq, wide, damping):
        if wide == 0:
            return (1 - damping) / (1 - eiaq * damping)
        return (1 - (eiaq * damping) ** wide) / (1 - eiaq * damping)
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the reflection and phase shift of the Film."""
        self._set_icq_(q, sinomega)
        r = super(self.__class__, self)._calc_reflection_(q, sinomega)
        r *= (self._calc_direction_(self._eiaq_ , self._wide_) * 
              self._calc_direction_(self._eibq_ , self._deep_) *
              self._calc_out_direction(self._eicq_ , self._high_, self._damping_)
             )
        return r      

    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the layer"""
        return [self._eiaq_ ** self._wide_,
                self._eibq_ ** self._deep_, 
                (self._eicq_ * self._damping_) ** self._high_
               ]
    
    def _calc_thickness_(self):
        """returns the thickness of the UnitCell."""
        if self._high_ == 0:
            return 1
        else:
            return self._high_ * self.crystal.c
    
    def _calc_depth_(self):
        """returns the depth of the UnitCell"""
        if self._deep_ == 0:
            return 1
        else:
            return self.crystal.b * self._deep_ 
    
    def _calc_width_(self):
        """returns the thickness of the UnitCell"""
        if self._wide_ == 0:
            return 1
        else:
            return self.crystal.a * self._wide_

########
# Film #
########

class Film(Cell):
    
    def __init__(self, cell, wide=0, deep=0, high=0, sigma = 0):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           layers:  number of layers in the film
        """
        self._cell_ = cell
        self._high_ = high
        self._deep_ = deep
        self._wide_ = wide
        self.eps = 0.0007*2*m.pi   
        self._Film__type_check__()
        self._sigma_ = sigma
        self._init_peaks_()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           cell -> Cell
           wide -> int
           deep -> int
           high -> int 
        """
        assert isinstance(self._cell_, (Cell)) , '_cell_ needs to be a Cell'
        assert isinstance(self._high_, (int)) , 'high needs to be a int'
        assert isinstance(self._deep_, (int)) , 'deep needs to be a int'
        assert isinstance(self._wide_, (int)) , 'wide needs to be a int'
        
    _Film__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self,  q, sinomega):
        """recalculates the values for i*(latticeparameters*q*n)"""
        factor = self.damping * self._cell_.c / sinomega
        self._iaq_ = -1J * self._cell_.a * q[0]
        self._ibq_ = -1J * self._cell_.b * q[1] 
        self._icq_ = -1J * self._cell_.c * q[2]
        self._eiaq_ = mc.exp(self._iaq_)
        self._eibq_ = mc.exp(self._ibq_)
        self._eicq_ = mc.exp(self._icq_)
        self._damping_ = m.exp(-factor)
    
    def _calc_direction_(self, eiaq, wide):
        #if wide == 0:
        if (eiaq.real > 1 - self.eps):
            return 1.0
            #return 0.0
        return self._calc_out_direction(eiaq, wide, 0.9999999999)
    
    def _calc_out_direction(self, eiaq, wide, damping):
        if wide == 0:
            return (1) / (1 - eiaq * damping)
        return (1 - (eiaq * damping) ** wide) / (1 - eiaq * damping)
    
    def _init_peaks_(self):
        self._q001_ = 2*m.pi/self._cell_.c
        self._I001_ = 1.0
        self._I002_ = 10.0
        
    def _roughness_(self, q):
        """Calculates a surface roughness gausian term"""
        if (self._sigma_ == 0.0): return 1.0
        #self._init_peaks_()
        return m.exp(-self._sigma_* (q[2] - self._q001_)**2)
    
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the reflection and phase shift of the Film."""
        r = self._cell_._calc_reflection_(q, sinomega)
        self._set_icq_(q, sinomega)
        r *= (self._calc_direction_(self._eiaq_ , self._wide_) * 
              self._calc_direction_(self._eibq_ , self._deep_) *
              self._calc_out_direction(self._eicq_ , self._high_, self._damping_) *
              self._roughness_(q)
             )
        return r

    def _calc_damping_(self):
        """returns the damping factor of one unitCell of the underlying crystal"""
        return self._cell_.damping

    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the Film"""
        return [self._eiaq_ ** self._wide_,
                self._eibq_ ** self._deep_, 
                (self._eicq_ * self._damping_) ** self._high_
               ]
    
    def _calc_thickness_(self):
        """returns the thickness of the UnitCell."""
        if self._high_ == 0:
            return 1
        else:
            return self._high_ * self._cell_.c
    
    def _calc_depth_(self):
        """returns the depth of the UnitCell"""
        if self._deep_ == 0:
            return 1
        else:
            return self._cell_.b * self._deep_ 
    
    def _calc_width_(self):
        """returns the thickness of the UnitCell"""
        if self._wide_ == 0:
            return 1
        else:
            return self._cell_.a * self._wide_        

#############
# Substrate #
#############

class Substrate(Film):

    def _calc_reflection_(self, q, sinomega):
        """calculates the reflection and phase shift of the Film."""
        r = self._cell_._calc_reflection_(q, sinomega)
        self._set_icq_(q, sinomega)
        r *= (self._calc_direction_(self._eiaq_ , self._wide_) * 
              self._calc_direction_(self._eibq_ , self._deep_) *
              (1) / ((self._q001_ - q[2])) *
              self._roughness_(q)
             )
        return r
##########
# Sample #
##########

class Sample(Cell):
        
    def __init__(self, films = None):
        """initializes the class data and checks data types
          
           Substrate: substrate of the Sample 
           Electrode: Electrode of the Sample
           Film: Film of the Sample
           Layers: any number of Layers to add to the sample
        """
        
        self._films_ = []
        if (films):
            for layer in films:
                self._films_.append(layer)
        self._Sample__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           _Layers_ -> list 
           _Layers_[i] -> Layer
        """
        assert isinstance(self._films_, list) , '_Layers_ needs to be an list'
        for i in range(len(self._films_)):
            assert (isinstance(self._films_[i], Cell)) , '_films_[' + str(i) + '] needs to be a Cell'
    
    _Sample__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def _calc_phase_shift_(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: scattering vector
           sinomega: sin of angle between surface and incedent beam
        """
        return [mc.exp(-1J * self.get_width() * q[0]), 
                mc.exp(-1J * self.get_depth() * q[1]), 
                mc.exp(-1J * self.get_thickness() * q[2]),
                ]
    
    def _phase_(self, film, q , sinomega):
        """returns the phaseshift of the film in the direction films are connected"""
        return 0
    
    def _calc_damping_(self):
        """returns the damping factor of one unitCell of the underlying cristal"""
        damping = 0
        for film in self._films_:
            damping += film.damping
        damping /= len(self._films_)
        return damping
    
    def _calc_reflection_(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           sinomega: sin of angle between surface and incedent beam
        """
        r = 0.0 + 0.0j        
        phase = 1.0+0.0j
        for layer in reversed(self._films_):
            r += layer.get_reflection(q, sinomega) * phase
            phase *= self._phase_(layer, q, sinomega)
        return r 
    
    def _calc_width_(self):
        """returns the width of the Sample."""
        thickness = 0.0
        for layer in self._films_:
            thickness += layer.get_width()
        thickness /= len(self._films_)
        return thickness
    
    def _calc_depth_(self):
        """returns the width of the Sample."""
        thickness = 0.0
        for layer in self._films_:
            thickness += layer.get_depth()
        thickness /= len(self._films_)
        return thickness
             
    
    def _calc_thickness_(self):
        """returns the thickness of the Sample."""
        thickness = 0.0
        for layer in self._films_:
            thickness += layer.get_thickness()
        thickness /= len(self._films_)
        return thickness
    
    def add_Layer(self, film):
        """adds the film to the top of the sample."""
        assert (isinstance(film, Cell)) , 'layer needs to be a Layer'
        self._films_.append(film)

#############
# MixSample #
#############

class MixSample(Sample):
    
    def __init__(self, films = None, probabileties = None):
        """initializes the class data and checks data types
          
           films: different films in the Sample
           probabileties: which part is made of which film
        """
        Sample.__init__(self,films)
        self._P_ = []
        if (probabileties):
            for p in probabileties:
                self._P_.append(p)
        self._MixSample__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           _P_ -> list 
           _P_[i] -> number
        """
        assert isinstance(self._P_, list) , '_P_ needs to be an list'
        for i in range(len(self._P_)):
            assert (isinstance(self._P_[i], _NUMBER_)) , '_P_[' + str(i) + '] needs to be a Number'
        self.__check_P__()
    
    _MixSample__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def __check_P__(self):
        """checks that the probabilities are still right"""
        eps = 0.0001
        assert (len(self._P_) == len(self._films_)) , 'each film needs a probability'
        if (self._P_):
            if (sum(self._P_) > 1.0+eps or sum(self._P_) < 1.0-eps):
                assert (sum(self._P_) == 1.0) , 'sum of all probabilities needs to be 1 but is ' + str(sum(self._P_))    
    
    def _calc_reflection_(self, q, sinomega):
        """returns the complex reflection value of the Sample.
             
           q: wave vector
           sinomega: sin of angle between surface and incedent beam
        """
        self.__check_P__()
        r = 0.0 + 0.0j
        i = 0        
        for layer in self._films_:
            r += layer.get_reflection(q, sinomega) * self._P_[i]
            i += 1
        return r
    
    def _calc_phase_shift_(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: scattering vector
           sinomega: sin of angle between surface and incedent beam
        """
        phase = [0.0,0.0,0.0]
        i = 0
        for layer in self._films_:
            for j in range(3):  
                phase[j] += layer.get_phase_shift(q, sinomega)[j] * self._P_[i]
            i += 1
        return phase
    
    def _calc_width_(self):
        """returns the width of the Sample."""
        thickness = 0.0
        i = 0
        for layer in self._films_:
            thickness += layer.get_width() * self._P_[i]
            i += 1
        return thickness
    
    def _calc_depth_(self):
        """returns the width of the Sample."""
        thickness = 0.0
        i = 0
        for layer in self._films_:
            thickness += layer.get_depth() * self._P_[i]
            i += 1
        return thickness
    
    def _calc_thickness_(self):
        """returns the thickness of the Sample."""
        t = 0
        i = 0
        for film in self._films_:
            t += film.get_thickness() * self._P_[i]
            i += 1
        return t  
    
##############
# SampleWide #
##############

class WideSample(Sample):
    
    def __init__(self, films = None):
        """initializes the class data and checks data types
          
           films: different films in the Sample
           probabileties: which part is made of which film
           widths: width of each film
        """
        super(self.__class__, self).__init__(films)
        
    def _phase_(self, film, q, sinomega):
        """returns the phaseshift of the film in the direction films are connected"""
        return film.get_phase_shift(q, sinomega)[0] 
    
    def _calc_width_(self):
        """returns the width of the Sample."""
        thickness = 0.0
        for layer in self._films_:
            thickness += layer.get_width()
        return thickness
    
###############
# ThickSample #
###############

class ThickSample(Sample):
    
    def __init__(self, films = None):
        """initializes the class data and checks data types
          
           films: different films in the Sample
           probabileties: which part is made of which film
           widths: width of each film
        """
        super(self.__class__, self).__init__(films)
        
    def _phase_(self, film, q, sinomega):
        """returns the phaseshift of the film in the direction films are connected"""
        return film.get_phase_shift(q, sinomega)[2]

    def _calc_thickness_(self):
        """returns the thickness of the Sample."""
        thickness = 0.0
        for layer in self._films_:
            thickness += layer.get_thickness()
        return thickness   
    
    def _calc_phase_shift_(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: scattering vector
           sinomega: sin of angle between surface and incedent beam
        """
        p3 = 1
        for layer in self._films_:
            p3 *= layer.get_phase_shift(q,sinomega)[2]
        return [mc.exp(-1J * self.get_width() * q[0]), 
                mc.exp(-1J * self.get_depth() * q[1]), 
                p3,
                ]
 
###################
# ComplexThikFilm #
###################

class ComplexThikFilm(MixSample):
    
    def __init__(self, cell, layers=1.0, width = 0.01, filmtype = Film):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           layers:  number of layers in the film
           width: width of the gausian distribution around layers
           filmType: the used filmtype to calculate the thinfilms
        """
        MixSample.__init__(self)
        self._cell_ = cell
        self._layers_ = layers
        self._width_ = width  
        self._filmType_ = filmtype
        self._ComplexThinFilm__type_check__()
        self._init_P_()        
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           layers -> number
           width -> number 
        """
        assert isinstance(self._layers_, _NUMBER_) , 'layers needs to be a Number'
        assert isinstance(self._width_, _NUMBER_) , 'width needs to be a Number'
        assert isinstance(self._cell_, Cell) , 'cell needs to be a Cell'
        #assert isinstance(self.filmType, (ThinFilm)) , 'filmType needs to be a ThinFilm'
        
    _ComplexThinFilm__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _init_P_(self):
        """initializes the probabilities for different Thicknesses"""
        self._P_ = []
        self._films_ = []
        minimum = int(self._layers_ - 3 * self._width_) #int always rounds down
        maximum = int(round(self._layers_ + 3 * self._width_ + 1.5)) #+0.5 to round up +1 for the range
        summe = 0
        for i in range(minimum, maximum):
            value = m.exp((i - self._layers_)**2/(-2 * self._width_**2))
            self._P_.append(value)
            self._films_.append(self._filmType_(self._cell_,0,0,i))
            summe += value
        for i in range(len(self._P_)):
            self._P_[i] /= summe
            
    def _init_P_perfect(self):
        """initializes the probabilities for different Thicknesses"""
        self._P_ = []
        self._films_ = []
        p2 = self._layers_ - int(self._layers_) #probabilety for second film
        p1 = 1-p2
        self._P_.append(p1)
        self._films_.append(Film(self._cell_,int(self._layers_)))
        self._P_.append(p2)
        self._films_.append(self._filmType_(self._cell_,int(self._layers_)+1))
 
    
##############
# BasicLayer #
##############

class BasicLayer(object):
        
    def __init__(self):
        """initializes the class data and checks data types"""
        pass
        
    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Basic Layer."""
        return 1
    
    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the Basic layer
        
            q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return [1,1,1]
        
    def _calc_thickness_(self):
        """returns the thickness of the Basic Layer"""
        return 1
    
    def _calc_N_(self):
        """returns the total number of layers of the Basic layer"""
        return 1
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Film.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return self._calc_reflection_(q, sinomega)
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the film.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        return self._calc_phase_shift_(q, sinomega)
    
    def get_thickness(self):
        """returns the thickness of the film."""
        return self._calc_thickness_()
    
    def get_N(self):
        """returns the total number of layers of the Basic layer"""
        return self._calc_N_()

###################
# BasicLayerCheck #
###################

class BasicLayerCheck(BasicLayer):
    
    def __init__(self):
        self._q_ = [0,0,0]
        self._sinomega_ = 0.0
        self._reflection_ = 0.0 + 0.0j
        self._phase_shift_ = [0.0 + 0.0j,0.0j,0.0j]
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Basic Layer.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        self._check_q_s_(q, sinomega)
        return self._reflection_
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the film.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        self._check_q_s_(q, sinomega)
        return self._phase_shift_
    
    def _check_q_s_(self, q, sinomega):
        """checks if q or sinomega are changed.
        
           q: scattering vector
           sinomega: sin(angle between incident beam and sample surface 
        """
        #assert isinstance(q,self._q_.__class__), 'q needs to be ' + str(self._q_.__class__)
        #assert isinstance(sinomega,self._sinomega_.__class__), 'sinomega needs to be ' + str(self._sinomega_.__class__)
        if ( (self._q_ != q) or 
             (self._sinomega_ != sinomega) ):
            self._q_ = q[:]
            self._sinomega_ = sinomega
            self._reflection_ = self._calc_reflection_(q, sinomega)
            self._phase_shift_ = self._calc_phase_shift_(q, sinomega)
            return True
        return False

##################
# BasicLayerSave #
##################

class BasicLayerSave(BasicLayer):
    
    def __init__(self):
        self._data_ = []
        self._reflection_ = []
        self._phase_shift_ = []
        self._counter_ = 0
        
    def _check_q_s_(self, q, sinomega):
        """checks if q or sinomega are changed.
        
           q: scattering vector
           sinomega: sin(angle between incident beam and sample surface""" 
        values = [q[:], sinomega]
        self._counter_ += 1
        try:
            if (values == self._data_[self._counter_]): #exception if counter bigger than list
                pass
            elif (values == self._data_[self._counter_-1]):
                self._counter_ -= 1
            else:                
                self._counter_ = self._data_.index(values) #exception if values not in data
        except:
            self._data_.append([q[:], sinomega])
            self._reflection_.append(self._calc_reflection_(q, sinomega))
            self._phase_shift_.append(self._calc_phase_shift_(q, sinomega))
            self._counter_ = (len(self._data_) - 2) #makes sure that every data set gets checked
            return self._counter_ + 1
        return self._counter_
        
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Basic Layer.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return self._reflection_[self._check_q_s_(q, sinomega)]
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the Basic Layer.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        return self._phase_shift_[self._check_q_s_(q, sinomega)][:]

#########
# Layer #
#########

class Layer(BasicLayer):
        
    def __init__(self, crystal):
        """initializes the class data and checks data types
          
           crystal CrystalStructure of the film 
        """
        BasicLayer.__init__(self)
        self.crystal= crystal
        self._icq_ = self._ibq_ = self._iaq_ = 0.0
        self._eiaq_ = self._eibq_ = self._eicq_ = 0.0
        self._damping_ = 0.0
        self.delta = 0.0
        
        self._Layer__type_check__()
    
    def __str__(self):
        """Returns the data in a structured form.
           
           returns its crystal structure 
        """
        return self.crystal.__str__()
    
    def __type_check__(self):
        """checks the class data to have the right types
        
           crystal -> Crystalstructure 
        """
        #assert isinstance(self.crystal, CrystalStructure) , 'crystal needs to be a SimpleCrystalstructure'
        
    _Layer__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self, q, sinomega):
        """recalculates the values for i*(latticeparameters*q)"""                   
        a = self.crystal.a
        b = self.crystal.b
        c = self.crystal.c
        factor = self.crystal.damping * c / sinomega
        self._iaq_ = -1J * a * q[0]
        self._ibq_ = -1J * b * q[1] 
        self._icq_ = -1J * c * q[2]
        self._eiaq_ = mc.exp(self._iaq_)
        self._eibq_ = mc.exp(self._ibq_)
        self._eicq_ = mc.exp(self._icq_)
        self._damping_ = mc.exp(-factor)
    
    def _calc_amplitude_(self, q):
        """returns calculates the scattering amplitude"""
        a = self.crystal.a
        b = self.crystal.b
        R = 2.1879e-15 #R= 2.1879e-15 is the classical electron radius
        amp1d = 4 * m.pi * R / a / b #/ m.sqrt(sum([x**2 for x in q]))
        return amp1d 
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Film."""
        self._set_icq_(q, sinomega)
        amp = self._calc_amplitude_(q)
        r = (amp * self.crystal.get_structure_factor(q, sinomega) / 
             ((1 - self._eiaq_ * self._damping_) * 
              (1 - self._eibq_ * self._damping_)
             )             
            )
        return r
    
    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the layer
        
            q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return [self._eiaq_, 
                self._eibq_, 
                mc.exp(self.delta * self._icq_)*self._eicq_*self._damping_
               ]
    
    def _calc_thickness_(self):
        """returns the thickness of the film."""
        return (1 + self.delta)*self.crystal.c

##############
# LayerCheck #
##############

class LayerCheck(Layer, BasicLayerCheck):
    
    def __init__(self, crystal):
        Layer.__init__(self, crystal)
        BasicLayerCheck.__init__(self)

#############
# LayerSave #
#############

class LayerSave(Layer, BasicLayerSave):
    
    def __init__(self, crystal):
        Layer.__init__(self, crystal)
        BasicLayerSave.__init__(self)          

#############
# ThickFilm #
#############
 
class ThickFilm (Layer):    
    
    def _set_icq_(self,  q, sinomega):
        """recalculates the values for i*(latticeparameters*q*n)"""
        super(self.__class__, self)._set_icq_(q, sinomega)
        self._damping2_ = self._damping_
        self._damping_ = 0.9999999
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Thick Film."""
        r = Layer._calc_reflection_(self, q, sinomega)
        r *= 1 / (1 - self._eicq_ * self._damping2_) * (1-self._damping_)**2
        #r = 2 * r / (1 + m.sqrt(1 + 4 * abs(r**2)))
        return r

##################
# ThickFilmCheck #
##################

class ThickFilmCheck (LayerCheck, ThickFilm):
    pass

#################
# ThickFilmSave #
#################

class ThickFilmSave (LayerSave, ThickFilm):
    pass

############
# ThinFilm #
############

class ThinFilm(Layer):

    def __init__(self, crystal, layers=0, breath=0, width=0):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           layers:  number of layers in the film
        """
        super(self.__class__, self).__init__(crystal)
        self._layers_ = layers
        self._breath_ = breath
        self._width_ = width
        self._ianq_ = self._ibnq_ = self._icnq_ = 0.0+0.0j
        self._eianq_ = self._eibnq_ = self._eicnq_ = 0.0+0.0j   
        self._ThinFilm__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           layers -> int 
        """
        assert isinstance(self._layers_, (int)) , 'layers needs to be a int'
        
    _ThinFilm__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self,  q, sinomega):
        """recalculates the values for i*(latticeparameters*q*n)"""
        super(ThinFilm, self)._set_icq_(q, sinomega)
        self._ianq_ =  self._breath_ * self._iaq_
        self._ibnq_ =  self._width_  * self._ibq_
        self._icnq_ =  self._layers_ * self._icq_
        self._eianq_ = self._eiaq_**self._breath_ if self._breath_ > 0 else 0 
        self._eibnq_ = self._eibq_**self._width_  if self._width_ > 0 else 0 
        self._eicnq_ = self._eicq_**self._layers_ if self._layers_ > 0 else 0 
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the reflection and phase shift of the Film."""
        r = Layer._calc_reflection_(self, q, sinomega)
        r *= ((1 - self._eianq_ * self._damping_ ** self._breath_) * 
              (1 - self._eibnq_ * self._damping_ ** self._width_)  *
              (1 - self._eicnq_ * self._damping_ ** self._layers_) / 
              (1 - self._eicq_ * self._damping_)
             )
        return r      

    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the layer"""
        return [self._eianq_ * self._damping_ ** self._breath_,
                self._eibnq_ * self._damping_ ** self._width_, 
                self._eicnq_ * self._damping_** self._layers_
               ]
    
    def _calc_thickness_(self):
        """returns the thickness of the film."""
        return (self._layers_ + self.delta)*self.crystal.c
    
    def _calc_N_(self):
        """returns the total number of layers of the Basic layer"""
        return self._layers_

#################
# ThinFilmCheck #
#################

class ThinFilmCheck(LayerCheck, ThinFilm):    
    
    def __init__(self, crystal, layers=1, breath=0, width=0):
        ThinFilm.__init__(self, crystal, layers, breath, width)
        LayerCheck.__init__(self, crystal)
        
################
# ThinFilmSave #
################

class ThinFilmSave(LayerSave, ThinFilm):    
    
    def __init__(self, crystal, layers=1, breath=0, width=0):
        ThinFilm.__init__(self, crystal, layers, breath, width)
        LayerSave.__init__(self, crystal)    

##############
# MultiLayer #
##############

class Multilayer(Sample):
    
    def __init__(self, film1 = None, film2 = None, film3 = None, films = None):
        Sample.__init__(self, substrate=film1, electrode=film2, film=film3, layers=films)
        
    def get_c(self):
        return self.get_thickness()/self.get_N()    
       
################
# SuperLattice #
################

class SuperLattice(Multilayer):
        
    def __init__(self, film1 = None, film2 = None, film3 = None, films = None, bilayers=1):
        
        """initializes the class data and checks data types
          
           film1: first film of the superlattice 
           film2: second film of the superlattice 
           bilayers: number of times both films get repeted
        """
        Sample.__init__(self, substrate=film1, electrode=film2, film=film3, layers=films)
        self.bilayers = bilayers
        self._SuperLattice__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           bilayers -> int 

        """
        assert isinstance(self.bilayers, (int)) , 'bilayers needs to be an int'
    
    _SuperLattice__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def _calc_phase_shift_(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        p = Multilayer._calc_phase_shift_(self, q, sinomega)
        return [p[0], p[1], p[2] ** self.bilayers]
    
    def _calc_reflection_(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        r = Multilayer._calc_reflection_(self, q, sinomega)
        p = Multilayer._calc_phase_shift_(self, q, sinomega)[2]
        reflection = 0.0 + 0.0j
        phase = 1.0
        for i in range(self.bilayers):
            reflection += (r * phase)
            phase *= p
        return reflection
    
    def _calc_thickness_(self):
        """returns the total thickness of the superlatice"""
        return self.get_multilayer_thickness() *self.bilayers
    
    def _calc_N_(self):
        """returns the total number of layers in the Superlattice"""
        return self.get_n() * self.bilayers
    
    def get_c(self):
        """returns the average c lattice parameter"""
        return self.get_multilayer_thickness() / self.get_n()
    
    def get_multilayer_thickness(self):
        """returns the Bilayer thickness of the superlatice"""
        return Multilayer._calc_thickness_(self)
    
    def get_n(self):
        """returns the number of layers per multilayer"""
        return Multilayer._calc_N_(self)
    

##############
# Experiment #
##############

class Experiment(object):
    
    def __init__(self, sample, base_struc, direct = 2.0E8, background = 1, wavelength = 1.5409E-10):
        self.sample = sample
        self.base_struc = base_struc
        self.direct = direct
        self.background = background
        self.wavelength = wavelength
        self._polarisation_ = True
        self._laue_ = True
    
    def gen_hkl(self, mini, maxi, step):
        hkl = mini
        while hkl < maxi:
            yield hkl
            hkl += step
            
    def gen_l(self, l_min = 0.8, l_max = 1.05, l_step = 0.0002, h = 0, k = 0, sinomegain = 0.0):
        l = l_min
        q = [0,0,0]
        q[0] = h * 2 * m.pi / self.base_struc.a
        q[1] = k * 2 * m.pi / self.base_struc.b
        factor = 2 * m.pi / self.base_struc.c
        _4pi = 1.0 / 4.0 / m.pi * self.wavelength
        while l < l_max:
            l += l_step
            q[2] = l * factor
            absq = m.sqrt(sum([x**2 for x in q]))
            if (sinomegain == 0.0):
                sinomega = absq  * _4pi
            else:
                sinomega = sinomegain
            yield [q,sinomega]
    
    def l_scan(self, l_min = 0.8, l_max = 1.05, l_step = 0.0002,  h = 0, k = 0, sinomegain = 0):

        return map(self._l_scan_, self.gen_l(l_min, l_max, l_step, h, k, sinomegain))
    
    def qz_scan(self, qz, sinomegain = 0):
        q = [0,0,0]
        data = []
        _4pi = 1.0 / 4.0 / m.pi * self.wavelength
        for point in qz:
            q[2] = point
            if (sinomegain == 0.0):
                sinomega = abs(point)  * _4pi
            else:
                sinomega = sinomegain
            data.append([q[:],sinomega])
        return map(self._l_scan_, data)
    
    def l_scan_multi(self, l_min = 0.8, l_max = 1.05, l_step = 0.0002,  h = 0, k = 0, sinomegain = 0):
             
        pool = multi.Pool()
        scans = pool.map(self._l_scan_, self.gen_l(l_min, l_max, l_step, h, k, sinomegain))
        pool.close()
        return scans
    
    def l_h_map(self, l_min = 0.8, l_max = 1.05, l_step = 0.0002,  h_min = -0.07, h_max = 0.07, h_step = 0.0007, k = 0, sinomegain = 0.0):
        
        data = []
        for h in self.gen_hkl(h_min, h_max, h_step):
            data.append(self.l_scan(l_min, l_max, l_step, h, k, sinomegain))
        return data
    
    def l_k_map(self, l_min = 0.8, l_max = 1.05, l_step = 0.0002,  k_min = -0.07, k_max = 0.07, k_step = 0.0007, h = 0, sinomegain = 0.0):
        
        data = []
        for k in self.gen_hkl(k_min, k_max, k_step):
            data.append(self.l_scan(l_min, l_max, l_step, h, k, sinomegain))
        return data
            
    def _l_scan_(self, inputs):
        q = inputs[0]
        sinomega = inputs[1]
        lin = abs(self.sample.get_reflection(q, sinomega))
        lin *= lin
        lin *= self.direct
        omeg = m.asin(sinomega) 
        if (self._laue_):
            lin *= 1.0/m.sin(2*omeg) #/sinomega
        if (self._polarisation_):
            lin *= (m.cos(2*omeg)**2) 
        lin += self.background
    
        return [q[0], q[1], q[2], lin] 

########
# Scan #
########

class Scan(object):
    
    def __init__(self, path = None, N1 = 2.0, N2 = 6.0, scans1 = 10.0, scans2 = 26.0, time1 = 90.0, time2 = 440.0, scan_time = 25.95, delay = 100.0 , data=None):
        self.N1 = N1
        self.N2 = N2
        self.scans1 = scans1
        self.scans2 = scans2
        self.time1 = time1
        self.time2 = time2
        self.scan_time = scan_time
        self.delay = delay
        if(data):
            self._data_ = data
        if (path):
            self._data_ = import_scan(path)
    
    def __getitem__(self, index):
        """implements the [] operator to get the data at index"""
        return self._data_[index]
    
    def __len__(self):
        return len(self._data_)
    
    def scan(self, index):
        """returns the scan index"""
        return self[index]
    
    def qz(self, index):
        """returns the qz valus of scan index"""
        return [float(a[0]*10**10) for a in self[index]]
    
    def int(self, index):
        """returns the intensety values of scan index"""
        return [float(a[1]) for a in self[index]]
    
    def _layer_(self, layer): 
        """returns the scan number after layer layers were grown on the sample"""
        Scans = self.scans1 + self.scans2
        N = self.N1 + self.N2
        scan = int(layer / N) * Scans
        
        temp = layer % N
        if (temp == 0):
            pass
        elif (temp > self.N1):
            scan += self.scans1
            scan += int( (self.delay + self.time2 * (temp - self.N1) / self.N2) / self.scan_time) + 1
        else:
            scan += int( (self.delay + self.time1 * temp / self.N1) / self.scan_time) + 1
        return scan
            
    def layer(self, layer): 
        """returns the scan after layer layers were grown on the sample"""
        return self[int(self._layer_(layer))]
    
    def qz_layer(self,layer):
        """returns the qz values after layer layers were grown on the sample"""
        return [float(a[0]*10**10) for a in self.layer(layer)]
    
    def int_layer(self,layer):
        """returns the intesety values after layer layers were grown on the sample""" 
        return [float(a[1]) for a in self.layer(layer)]
    
    
##################
# Help Functions #
##################

def l_scan(sample, base_struc, l_min = 0.8, l_max = 1.05, l_step = 0.0002,  h = 0, k = 0, sinomegain = 0, direct = 2.0E8, background = 1, wavelength = 1.5409E-10):

    q = [0.0,0.0,0.0]
    steps_l = (l_max-l_min)/l_step
    q[0] = h*2*m.pi/base_struc.a
    q[1] = k*2*m.pi/base_struc.b
    scan = []
    for i in range(int(steps_l)):
        l=l_min+i*l_step
        q[2] = l*2*m.pi/base_struc.c
        absq = m.sqrt(sum([x**2 for x in q]))
        if (sinomegain == 0):
            sinomega = absq * wavelength / 4 / m.pi
        else:
            sinomega = sinomegain
        lin = abs(sample.get_reflection(q, sinomega))
        lin *= lin
        lin *= direct
        lin += background
        #lin *= 1/m.sin(2*m.asin(sinomega))*(1+m.cos(2*m.asin(sinomega))**2)/2/sinomega
        scan.append([q[0], q[1], q[2], lin])
    
    return scan

def cen(l, inte):
    l_cen = 0
    for x,y in zip(l,inte):
        l_cen += x*y
    l_cen /= sum(inte)
    return l_cen

def get_values(data_in,size=[4,200,550,-1],width=10,N=8.0):
    data = []
    i = 0
    for scan in data_in:
        inte = [y[3] for y in scan]
        l = [y[2] for y in scan]
        int2 = max(inte[size[2]:size[3]])
        int1 = max(inte[size[1]:size[2]])
        int0 = max(inte[size[0]:size[1]])
        l0_in = inte.index(int0)
        l1_in = inte.index(int1)
        l2_in = inte.index(int2)
        if (l0_in-width < size[0]):
            l0 = cen(l[size[0]:size[0] + 2 * width + 1],inte[size[0]:size[0] + 2 * width+ 1])
        else:
            l0 = cen(l[l0_in-width:l0_in + width + 1],inte[l0_in-width:l0_in + width + 1])
        l1 = cen(l[l1_in-width:l1_in + width + 1],inte[l1_in-width:l1_in + width + 1])
        l2 = cen(l[l2_in-width:l2_in + width + 1],inte[l2_in-width:l2_in + width + 1])
        c0 = 2*m.pi/l0
        c12 = 2*m.pi/(2*l1-l2)
        L01 = 2*m.pi/(l0-l1)
        L12 = 2*m.pi/(l1-l2)
        N011 = L01/c0
        N121 = L12/c0
        N012 = L01/c12
        N122 = L12/c12
        data.append([i,i/N,l0,l1,l2,c0,c12,L01,L12,N011,N012,N121,N122,int0,int1,int2])
        i += 1
    return data

def lcm(a,b): 
    return abs(a * b) / fractions.gcd(a,b) if a and b else 0 

def get_multi(x):
    S = str(x).split(".")
    if (len(S) == 1):
        return 1
    L = 10**len(S[1])
    div = alldividers(L)
    for i in div:
        if(i*x == float(int(i*x))):
            return i
                      
def alldividers(number):
    multi = 10.0
    try: #only changes the multiplicator if number is a decimal number
        temp = str(number).split(".")
        multi = 10.0**len(temp[1])
    except:
        pass        
    result = []
    for i in range (1, number + 1):
        if number * multi / i % multi == 0:
            result.append(i)
    return result

def create_layer(struc1, struc2, amount1):
    L = get_multi(amount1) 
    a = max(struc1.a,struc2.a)
    b = max(struc1.b,struc2.b)
    c = max(struc1.c,struc2.c)
    crystal = CrystalStructure([a*L,b,c])
    for j in range(L):
        if (j < amount1*L):
            struc = struc1                
        else:
            struc = struc2 
        for atom in struc.atoms:
            addatom = Atom(atom.form, atom.pos[:] )
            addatom.pos[0] += j*a
            crystal.add_atom(addatom)
    return crystal

def create_structure(a_in, c_in, a11 ,a12 ,a21 ,a22, l1, l2, bl, size = 0.0, zoffset1 = 0.0, zoffset2 = 0.0):
    random.seed()
    L1 = get_multi(l1)
    L2 = get_multi(l2)
    L = lcm(L1,L2)
    layers = bl*(l1+l2)
    if (float(int(layers)) == layers):
        layers = int(layers)
    else:
        layers = int(layers)+1
    (a, b, c) = (a_in*L*10**(-10), a_in*10**(-10) ,c_in*layers*10**(-10) )
    crystal = CrystalStructureCheck()
    crystal.a = a
    crystal.b = b
    crystal.c = c
    
    atom_kind = 1
    l = l1
    offset = 0
    zoffset = zoffset1
    for i in range(layers):
        for j in range(L):
            #check that the right atom is printed
            if (i*L+j - offset >= l*L):
                if (atom_kind == 1):
                    l = l2
                    atom_kind = 2
                    zoffset += zoffset2
                else:
                    l = l1
                    atom_kind = 1
                offset = i*L+j
                zoffset += zoffset1
            #create atoms
            #corner atom
            if (atom_kind == 1):
                atomtype = a11
            else:
                atomtype = a21
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/" + atomtype + ".at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.0 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.0 + random.randrange(-1,1) * size) 
            pos[2] =(i + 1.0 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #center atom
            if (atom_kind == 1):
                atomtype = a12
            else:
                atomtype = a22
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/" + atomtype + ".at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 1
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 1.0 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 2
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.0 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 3
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.0 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
    return crystal

def showup(l1, l2):
    L1 = get_multi(l1)
    L2 = get_multi(l2)
    L = lcm(L1,L2)
    bl = get_multi((l1+l2))
    layers = float(bl)*(l1+l2)
    
    if (float(int(layers)) == layers):
        layers = int(layers)  
    
    atom = 1
    l = l1
    offset = 0
    result = ''
    for i in range(layers):
        for j in range(L):
            #check that the right atom is printed
            if (i*L+j - offset >= l*L):
                if (atom == 1):
                    l = l2
                    atom = 2
                else:
                    l = l1
                    atom = 1
                offset = i*L+j
            if (atom == 1):
                result += 'B'
            else:
                result += 'S'
        result += '\n'
    each = result.split('\n')
    temp = []
    for i in range(len(each)-1, -1 ,-1):
        temp.append(each[i])
    #for line in temp:
        #print line
        
#######################
# Load Data Functions #
#######################

def load_sim(path):
    data = []
    f = open(path, 'r')
    line=f.readline()
    line=f.readline()
    while (line):
        words = line.split(' ')
        if words:
            if (len(words) > 1):
                data.append([words[2],words[3]])
        line=f.readline()
    f.close()
    return data

def load_data(path, row1 = 0, row2 = 1, separator=' '):
    data = []
    f = open(path, 'r')
    line=f.readline()
    line=f.readline()
    while (line):
        words = line.split(separator)
        if words:
            if (len(words) > 1):
                data.append([words[row1],words[row2]])
        line=f.readline()
    f.close()
    return data

def load_scan(path, scan_number):
    data = []
    f = open(path, 'r')
    line = f.readline()
    line = f.readline()
    while (line):
        words = line.split(',')
        if words:
            if (len(words) > 1):
                if (words[0] == str(scan_number)):
                    data.append([float(words[1]),float(words[3])])
        line=f.readline()
    f.close()
    return data

def import_scan(path, row1 = 1, row2 = 3, separator=','):
    data = []
    data.append([])
    f = open(path, 'r')
    line = f.readline()
    line = f.readline()
    offset = 0
    while (line):
        words = line.split(separator)
        if words:
            if (len(words) > 1):
                try:
                    data[int(words[0])-offset].append([float(words[row1]),float(words[row2]),float(words[0])])  
                except:
                    data.append([])
                    if (int(words[0])-offset >= len(data)):
                        offset = int(words[0])-len(data)+1
                    data[int(words[0])-offset].append([float(words[row1]),float(words[row2]),float(words[0])])
        line=f.readline()
    f.close()
    return data

def get_c2(c,c1, n1 = 6, n2 = 2):
    return (c * (n2 + n1) - n1 * c1) / n2 

