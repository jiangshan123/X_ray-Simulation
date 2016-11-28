#include <math.h>
#include "3D_math.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <complex>

struct Atom;
struct StructureFactor;
struct Perovskite;
struct CrystalStructure;
class Film;
class ThikFilm;
class Sample;

//////////
// Atom //
//////////

struct Atom {
	Vec3D pos;
	StructureFactor* struc;
	Atom* next;
	Atom(double x, double y, double z, StructureFactor* struc);
	Atom(Vec3D position, StructureFactor* struc);
	Atom();
	void Free();
};

//////////////////////
// structure factor //
//////////////////////

struct StructureFactor {
	Vec3D q;
	double a[5];
	double b[4];
	std::string atom;
	double value;
	double get_value(Vec3D q);
	int load_atom(std::string path);
	StructureFactor(std::string atom, double a1, double a2, double a3, double a4, double a5, double b1, double b2, double b3, double b4);
	StructureFactor(std::string path);
	StructureFactor();
};

std::ostream& operator<<( std::ostream& os, const StructureFactor& S );

////////////////
// Perovskite //
////////////////

struct Perovskite {
	Vec3D lattice_parameters;
	Vec3D q;
	StructureFactor* A_structurfactor;
	StructureFactor* B_structurfactor;
	StructureFactor* O_structurfactor;
	std::complex<double> form_factor;
	std::complex<double> get_form_factor(Vec3D q);
	Perovskite(Vec3D lattice_parameters, StructureFactor* A_structurfactor, StructureFactor* B_structurfactor, StructureFactor* O_structurfactor);
	Perovskite();
};

std::ostream& operator<<( std::ostream& os, const Perovskite& P );

//////////////////////
// crystalstructure //
//////////////////////

struct CrystalStructure{
	Atom* atoms;
	Vec3D lattice_parameters;
	Vec3D q;
	std::complex<double> form_factor;
	std::complex<double> get_form_factor(Vec3D q);
	CrystalStructure(Vec3D lattice_parameters, Atom* atoms);
	CrystalStructure(std::string path);
	CrystalStructure();
	int load_structure(std::string path);
};

std::ostream& operator<<( std::ostream& os, const Perovskite& P );

//////////
// Film //
//////////

class Film {

public:
	//Constructors
	Film();
	Film(CrystalStructure* crystal, double damping);
	//Destructors
	~Film();
	//Functions
	std::complex<double> get_reflection(Vec3D q, double wavelength);
	std::complex<double> get_reflection_sum(Vec3D q, double wavelength, int n);
	void set_icq(Vec3D q, double wavelength);
	void calc_amplitude(Vec3D q);
protected:
	CrystalStructure* crystal;
	double damping, wavelength;
	Vec3D q;
	std::complex<double> iaq, ibq, icq;
	std::complex<double> reflection;
	std::complex<double> amplitude;
private:
	void Init();
};

//////////////
// ThinFilm //
//////////////

class ThinFilm : Film {
public:
	//Constructors
	ThinFilm();
	ThinFilm(CrystalStructure* crystal, double layers, double damping);
	//Destructors
	~ThinFilm();
	//Functions
	std::complex<double> get_reflection(Vec3D q, double wavelength);
	std::complex<double> get_phase_shift(Vec3D q, double wavelength);
	void set_icnq(Vec3D q, double wavelength);
private:
	double layers;
	std::complex<double> icnq;
	std::complex<double> phase_shift;
	void Init();
	bool check_q_l(Vec3D q, double wavelength);
};

//////////////
// ThikFilm //
//////////////

class ThikFilm : Film {

public:
	//Constructors
	ThikFilm();
	ThikFilm(Perovskite* crystal, double damping);
	//Destructors
	~ThikFilm();
	//Functions
	std::complex<double> get_reflection(Vec3D q, double wavelength);
private:

};

////////////
// Sample //
////////////

class Sample {

private:
	StructureFactor* structurfactors;
	void initstructurfactors();

public:
	//Constructors
	Sample();
	//Destructors
	~Sample();
};
