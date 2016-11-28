#include "main.h"
#include <ctime>


int main(int argc, char* argv[])

{
	clock_t t;
	t = clock();
	std::string experiment_path = std::string(argv[1]);

	//scan constants
	double h_min = 0.0;
	double h_max = 0.0005;
	double h_step = 0.0005;

	double k_min = -0.0;
	double k_max = 0.0005;
	double k_step = 0.0005;

	double l_min = 0.8;
	double l_max = 1.05;
	double l_step = 0.0005;

	double lambda = 1.5409E-10;
	double direct = 1.5E8;
	double background = 0;
	//film constants
	double damping = 5.98E+4;
	double layers = 50;

	int unitcells = 1;
	std::string film_path;
	//returns 0 if successful
	//returns 1 if no path was given
	//returns 2 if no valid path was given
	int check = 3;
	std::ifstream ifile;
	if (experiment_path.c_str()) {
		ifile.open(experiment_path.c_str());
		if (ifile.is_open()) {
			std::string buffer;
			double value = 0;
			while (ifile.good() ) {
				ifile >> buffer;
				if (buffer == "h_min") {
					ifile >> value;
					h_min = value;
				}
				if (buffer == "h_max") {
					ifile >> value;
					h_max = value;
				}
				if (buffer == "h_step") {
					ifile >> value;
					h_step = value;
				}
				if (buffer == "k_min") {
					ifile >> value;
					k_min = value;
				}
				if (buffer == "k_max") {
					ifile >> value;
					k_max = value;
				}
				if (buffer == "k_step") {
					ifile >> value;
					k_step = value;
				}
				if (buffer == "l_min") {
					ifile >> value;
					l_min = value;
				}
				if (buffer == "l_max") {
					ifile >> value;
					l_max = value;
				}
				if (buffer == "l_step") {
					ifile >> value;
					l_step = value;
				}
				if (buffer == "lambda") {
					ifile >> value;
					lambda = value;
				}
				if (buffer == "direct") {
					ifile >> value;
					direct = value;
				}
				if (buffer == "background") {
					ifile >> value;
					background = value;
				}
				if (buffer == "layers") {
					ifile >> value;
					layers = value;
				}
				if (buffer == "damping") {
					ifile >> value;
					damping = value;
				}
				if (buffer == "film") {
					ifile >> buffer;
					film_path = buffer;
				}
				if (buffer == "unitcells") {
					ifile >> value;
					unitcells = (int) value;
				}
			}
			ifile.close();
			check = 0;
		} else {
			check = 2;
		}
	} else {
		check = 1;
	}

	if (check != 0) {
		std::cout << "error loading experiment error code: " << check << " from path" << experiment_path << std::endl;
	} else {

	double steps_h = (h_max-h_min)/h_step;
	double steps_k = (k_max-k_min)/k_step;
	double steps_l = (l_max-l_min)/l_step;
	//Load data
	/*StructureFactor* Sr = new StructureFactor("./atoms/Sr.at");
	StructureFactor* Ru = new StructureFactor("./atoms/Ru.at");
	StructureFactor* Pb = new StructureFactor("./atoms/Pb.at");
	StructureFactor* O = new StructureFactor("./atoms/O.at");
	StructureFactor* Ti = new StructureFactor("./atoms/Ti.at");*/

	//create crystals
	CrystalStructure* PTO = new CrystalStructure("./Perovskites/PTOs.str");
	CrystalStructure* PTO2 = new CrystalStructure(film_path);
	CrystalStructure* STO = new CrystalStructure("./Perovskites/STO.str");

	//create Films
	Film* substrate = new Film(STO,damping);
	ThinFilm* film = new ThinFilm(PTO,unitcells,damping);
	ThinFilm* film2 = new ThinFilm(PTO2,unitcells,damping);

	//Output
	Vec3D q;
	std::cout << "#h k l value log(value)" << std::endl;
	for (int i = 0; i<steps_h; i++) {
		double h=h_min+i*h_step;
		for (int j = 0; j<steps_k; j++) {
			double k=k_min+j*k_step;
			for (int m = 0; m<steps_l; m++) {
				double l=l_min+m*l_step;
				//q[1] = h*2*M_PI/STO->lattice_parameters[1];
				q[2] = l*2*M_PI/STO->lattice_parameters[2];
				//double lin = abs(film2->get_reflection(q, lambda)+substrate->get_reflection(q, lambda)*film2->get_phase_shift(q, lambda));
				//double lin = abs(film->get_reflection(q, lambda));
				double lin = abs(substrate->get_reflection(q, lambda));
				//double lin = abs(PTO2->get_form_factor(q));
				lin *= lin;
				lin *= direct;
				lin += background;
				//std::cout << h << " " << k << " " << l << " " << lin << " " << log10(lin) << std::endl;
			}
			//std::cout << std::endl;
		}
	}
	t = clock() - t;
	printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}}
