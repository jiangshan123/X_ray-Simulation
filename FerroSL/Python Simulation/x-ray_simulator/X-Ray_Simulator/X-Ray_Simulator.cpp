#include "X-Ray_Simulator.h"

//////////
// Atom //
//////////

//Intilization with 0 values
Atom::Atom() {
	next = NULL;
	struc = NULL;
}

//Intilization with given values
Atom::Atom(double x, double y, double z, StructureFactor* struc) {
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	next = NULL;
	this->struc = struc;
}

//Intilization with given values
Atom::Atom(Vec3D position, StructureFactor* struc) {
	this->pos = position;
	this->struc = struc;
	this->next = NULL;
}

//Destroys all atoms in the list
void Atom::Free() {
	if(this->next) this->next->Free();
	delete this;
}
//////////////////////
// structure factor //
//////////////////////

//Intilization with 0 values
StructureFactor::StructureFactor() {
	this->atom = "";
	this->a[0] = 0;
	this->a[1] = 0;
	this->a[2] = 0;
	this->a[3] = 0;
	this->a[4] = 0;
	this->b[0] = 0;
	this->b[1] = 0;
	this->b[2] = 0;
	this->b[3] = 0;
	this->value = 0;
}

//Intilization with given values
StructureFactor::StructureFactor(std::string atom, double a1, double a2, double a3, double a4, double a5, double b1, double b2, double b3, double b4) {
	this->atom = atom;
	this->a[0] = a1;
	this->a[1] = a2;
	this->a[2] = a3;
	this->a[3] = a4;
	this->a[4] = a5;
	this->b[0] = b1;
	this->b[1] = b2;
	this->b[2] = b3;
	this->b[3] = b4;
	this->value = 0;
}

//Intilization with values from file at path
StructureFactor::StructureFactor(std::string Path) {
	this->atom = "";
	this->a[0] = 0;
	this->a[1] = 0;
	this->a[2] = 0;
	this->a[3] = 0;
	this->a[4] = 0;
	this->b[0] = 0;
	this->b[1] = 0;
	this->b[2] = 0;
	this->b[3] = 0;
	this->value = 0;
	int check = this->load_atom(Path);
	if (check != 0) std::cout << "error loading Atom error code: " << check << " from path" << Path << std::endl;
}

//calculate the Structure Factor for a given q and returns it
double StructureFactor::get_value(Vec3D q) {
	if (this->q != q) {
		double q2 = 16*M_PI*M_PI/(q*q);
		this->value = a[0]*exp(b[0]*q2)+a[1]*exp(b[1]*q2)+a[2]*exp(b[2]*q2)+a[3]*exp(b[3]*q2)+a[4];
		this->q = q;
	}
	return this->value;
}

//loads constants for the structurefactor from a file
//returns 0 if succesfull
//returns 1 if no path was given
//returns 2 if no valid path was given
int StructureFactor::load_atom(std::string Path) {
	std::ifstream ifile;
	if (Path.c_str()) {
		ifile.open(Path.c_str());
		if (ifile.is_open()) {
			std::string buffer;
			double value = 0;
			while (ifile.good() ) {
				ifile >> buffer;
				if (buffer == "atom") {
					ifile >> buffer;
					this->atom = buffer;
				}
				if (buffer == "a1") {
					ifile >> value;
					this->a[0] = value;
				}
				if (buffer == "a2") {
					ifile >> value;
					this->a[1] = value;
				}
				if (buffer == "a3") {
					ifile >> value;
					this->a[2] = value;
				}
				if (buffer == "a4") {
					ifile >> value;
					this->a[3] = value;
				}
				if (buffer == "a5") {
					ifile >> value;
					this->a[4] = value;
				}
				if (buffer == "b1") {
					ifile >> value;
					this->b[0] = value;
				}
				if (buffer == "b2") {
					ifile >> value;
					this->b[1] = value;
				}
				if (buffer == "b3") {
					ifile >> value;
					this->b[2] = value;
				}
				if (buffer == "b4") {
					ifile >> value;
					this->b[3] = value;
				}
			}
			ifile.close();
			return 0;
		}
		return 2;
	}
	return 1;
}

//Ofstream operator for StructureFactor
std::ostream& operator<<( std::ostream& os, const StructureFactor& Struc ) {
	os << "element:" << Struc.atom << std::endl;
	for (int i=0; i<5; i++) {
		os << 'a' << i+1 << " = " << Struc.a[i] << "; ";
	}
	os << std::endl;
	for (int i=0; i<4; i++) {
			os << 'b' << i+1 << " = " << Struc.b[i] << "; ";
		}
	os << std::endl;
	return os;
}

//////////////////////
// CrystalStructure //
//////////////////////

//loads constants for the Crystal Structure from a file
//returns 0 if succesfull
//returns 1 if no path was given
//returns 2 if no valid path was given
int CrystalStructure::load_structure(std::string Path) {
	if (this->atoms) {
		this->atoms->Free();
		this->atoms = NULL;
	}
	std::ifstream ifile;
	if (Path.c_str()) {
		ifile.open(Path.c_str());
		if (ifile.is_open()) {
			std::string buffer;
			Atom* before = NULL;
			Atom* current = NULL;
			double value = 0;
			while (ifile.good() ) {
				ifile >> buffer;
				if (buffer == "a") {
					ifile >> value;
					this->lattice_parameters[0] = value;
				}
				if (buffer == "b") {
					ifile >> value;
					this->lattice_parameters[1] = value;
				}
				if (buffer == "c") {
					ifile >> value;
					this->lattice_parameters[2] = value;
				}
				if (buffer == "element") {
					if (this->atoms) {
						before = current;
						current = new Atom();
						before->next = current;
					} else {
						this->atoms = new Atom();
						current = this->atoms;
					}
					ifile >> buffer;
					current->struc = new StructureFactor("./atoms/"+buffer+".at");
				}
				if (buffer == "x") {
					ifile >> value;
					current->pos[0] = value*this->lattice_parameters[0];
				}
				if (buffer == "y") {
					ifile >> value;
					current->pos[1] = value*this->lattice_parameters[1];
				}
				if (buffer == "z") {
					ifile >> value;
					current->pos[2] = value*this->lattice_parameters[2];
				}
			}
			ifile.close();
			return 0;
		}
		return 2;
	}
	return 1;
}

//Intilization with 0 values
CrystalStructure::CrystalStructure() {
	this->atoms = NULL;
}

//Intilization with given values
CrystalStructure::CrystalStructure(Vec3D lattice_parameters, Atom* atoms) {
	this->atoms = atoms;
	this->lattice_parameters = lattice_parameters;
}

//Intilization with given values
CrystalStructure::CrystalStructure(std::string Path) {
	this->atoms = NULL;
	int check = load_structure(Path);
	if (check != 0) std::cout << "error loading Crystal Structure error code: " << check << " from path" << Path << std::endl;
}

//Calculates the formfactor of the Crystal Structure for the given q value
std::complex<double> CrystalStructure::get_form_factor(Vec3D q) {
	if (atoms) {
		if (q != this->q) {
			std::complex<double> i;
			i.imag(-1);
			Atom * current = this->atoms;
			this->form_factor.real(0);
			this->form_factor.imag(0);
			while (current) {
				this->form_factor += current->struc->get_value(q)*exp(i*(current->pos*q));
				current = current->next;
			}
		this->q =q;
		}
		return this->form_factor;
	}
	return 0;
}

////////////////
// Perovskite //
////////////////

//Intilization with given values
Perovskite::Perovskite(Vec3D lattice_parameters, StructureFactor* A_structurfactor, StructureFactor* B_structurfactor, StructureFactor* O_structurfactor) {
	this->lattice_parameters = lattice_parameters;
	this->A_structurfactor = A_structurfactor;
	this->B_structurfactor = B_structurfactor;
	this->O_structurfactor = O_structurfactor;
}

//Intilization with 0 values
Perovskite::Perovskite() {
	this->A_structurfactor = NULL;
	this->B_structurfactor = NULL;
	this->O_structurfactor = NULL;
}

//Calculates the formfactor of the perovskite for the given q value
std::complex<double> Perovskite::get_form_factor(Vec3D q) {
	if (this->A_structurfactor != NULL and this->B_structurfactor != NULL and this->O_structurfactor != NULL) {
		if (q != this->q) {
			std::complex<double> icq, iaq, ibq;
			//create list of atom position in unit cell and loop over each atom with
			//icq = I*q*atom position
			//form_factor and form_factor = structurfactor_atom*exp(-icq)
			iaq.imag(q[0]*this->lattice_parameters[0]);
			ibq.imag(q[1]*this->lattice_parameters[1]);
			icq.imag(q[2]*this->lattice_parameters[2]);
			this->form_factor = A_structurfactor->get_value(q)+O_structurfactor->get_value(q)*exp(-(iaq+ibq)*0.5) + B_structurfactor->get_value(q)*exp(-(iaq+ibq)*0.5-icq*0.5) +  O_structurfactor->get_value(q)*exp(-(iaq)*0.5-icq*0.5)+ O_structurfactor->get_value(q)*exp(-(ibq)*0.5-icq*0.5);
			//std::cout << "A: " << A_structurfactor->get_value(q) << " B: " << B_structurfactor->get_value(q) << " O: "  << O_structurfactor->get_value(q) <<  this->form_factor << std::endl;
			this->q =q;
		}
		return this->form_factor;
	}
	return 0;
}

//Ofstream operator for Perovskite
std::ostream& operator<<( std::ostream& os, const Perovskite& P ) {
	os << P.A_structurfactor->atom+P.B_structurfactor->atom+P.O_structurfactor->atom+"3"  << std::endl;
	return os;
}

//////////
// Film //
//////////

//Intilization with 0 values
Film::Film() {
	this->crystal = NULL;
	this->damping = 0;
	Init();
}

//Intilization with given values
Film::Film(CrystalStructure* crystal, double damping) {
	this->crystal = crystal;
	this->damping = damping;
	Init();
}

void Film::Init()  {
	this->wavelength = 0;
	this->icq = ibq = iaq = 0;
	this->reflection = 0;
	this->amplitude = 0;
}

void Film::set_icq(Vec3D q, double wavelength) {
	//damping???
	iaq.imag(-1*(this->crystal->lattice_parameters[0]*q[0]));
	iaq.real(-this->crystal->lattice_parameters[0]*4*M_PI/wavelength*this->damping/q.module());
	ibq.imag(-1*(this->crystal->lattice_parameters[1]*q[1]));
	ibq.real(-this->crystal->lattice_parameters[1]*4*M_PI/wavelength*this->damping/q.module());
	icq.imag(-1*(this->crystal->lattice_parameters[2]*q[2]));
	icq.real(-this->crystal->lattice_parameters[2]*4*M_PI/wavelength*this->damping/q.module());
}

std::complex<double> Film::get_reflection(Vec3D q, double wavelength) {
	if ( (this->q != q) or (this->wavelength != wavelength) ) {
			this->q =q;
			this->wavelength = wavelength;
			set_icq(q, wavelength);
			double one = 1;
			double two = 2;
			calc_amplitude(q);
			reflection = amplitude*crystal->get_form_factor(q)/((one-exp(iaq))*(one-exp(ibq))*(one-exp(icq)));
			reflection = two*reflection/(one+sqrt(one+two*two*(reflection.real()*reflection.real()+reflection.imag()*reflection.imag())));
		}
	return reflection;

}

std::complex<double> Film::get_reflection_sum(Vec3D q, double wavelength, int n) {
	if ( (this->q != q) or (this->wavelength != wavelength) ) {
			this->q =q;
			this->wavelength = wavelength;
			set_icq(q, wavelength);
			//calc_amplitude(q);
			reflection.real(0);
			reflection.imag(0);
			for (int i = 0; i < n; i++) {
				//for (int j = 0; j < n; j++) {
					//for (int k = 0; k < n; k++) {
						reflection += exp(icq*(double)i);//+ibq*(double)j+icq*(double)k);
					//}
				//}
			}
			reflection /= n;
			//reflection = amplitude*crystal->get_form_factor(q)/((one-exp(iaq))*(one-exp(ibq))*(one-exp(icq)));
			//reflection = two*reflection/(one+sqrt(one+two*two*(reflection.real()*reflection.real()+reflection.imag()*reflection.imag())));
		}
	return reflection;
}

//calculates the amplitude for the reflection
void Film::calc_amplitude(Vec3D q) {
	//R= 2.1879e-15 is the classical electron radius
	//where is a^2 comming from?
	double amplitude1d = 4*M_PI*2.1879e-15/q.module()/crystal->lattice_parameters[0]/crystal->lattice_parameters[0];
	amplitude.imag(amplitude1d*amplitude1d*amplitude1d);
}

//////////////
// ThinFilm //
//////////////

//Intilization with 0 values
ThinFilm::ThinFilm() : Film() {
	this->layers = 0;
	Init();
}

//Intilization with given values
ThinFilm::ThinFilm(CrystalStructure* crystal, double layers, double damping):Film(crystal, damping) {
	this->layers = layers;
	Init();
}

void ThinFilm::Init() {
	this->icnq = 0;
	this->phase_shift = 0;
}

//calculates the film diffraction
std::complex<double> ThinFilm::get_reflection(Vec3D q, double wavelength){
	check_q_l(q, wavelength);
	//std::cout << this->layers << " " << exp(nicq) << " " << exp(icq) << std::endl;
	return reflection;
}

//calculates the phase of the film
std::complex<double> ThinFilm::get_phase_shift(Vec3D q, double wavelength) {
	check_q_l(q, wavelength);
	return phase_shift;
}

//sets the icnq value of the film
void ThinFilm::set_icnq(Vec3D q, double wavelength) {
	//damping???
	icnq.imag(-this->layers*(this->crystal->lattice_parameters[2]*q[2]));
	icnq.real(-this->layers*this->crystal->lattice_parameters[2]*4*M_PI/wavelength*this->damping/q.module());
}

//checks if q or lambda are changed and recalcultes intensety and phase
bool ThinFilm::check_q_l(Vec3D q, double wavelength) {
	if ( (this->q != q) or (this->wavelength != wavelength) ) {
		this->q =q;
		this->wavelength = wavelength;
		set_icq(q, wavelength);
		set_icnq(q, wavelength);
		phase_shift = exp(icnq);
		calc_amplitude(q);
		double one = 1;
		reflection = amplitude*crystal->get_form_factor(q)*(one-phase_shift)/((one-exp(iaq))*(one-exp(ibq))*(one-exp(icq)));
		return false;
	}
	return true;
}
