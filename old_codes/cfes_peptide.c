#include "Start_C.h"

int main(int argc, char **argv)
{
//LAMMPS initialization................................
	std::string directory = argv[1];
	Start(argc, argv);

//Files.................................................
	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps"));
	lammps_file(LAMMPS_HANDLE, lammps_input_file_address.c_str());

	std::string input_file_address = getFileAddress(directory, std::string("input.in"));
	std::ifstream input_file(input_file_address.c_str());

	std::string hist_file_address = getFileAddress(directory, std::string("hist.dat"));
	std::ifstream hist_file(hist_file_address.c_str());

	std::string xcoord_file_address = getFileAddress(directory, std::string("xcoord.dat"));
	std::ifstream xcoord_file(xcoord_file_address.c_str());

	std::string ycoord_file_address = getFileAddress(directory, std::string("ycoord.dat"));
	std::ifstream ycoord_file(ycoord_file_address.c_str());

	std::string zcoord_file_address = getFileAddress(directory, std::string("zcoord.dat"));
	std::ifstream zcoord_file(zcoord_file_address.c_str());

	std::string force_x_file_address = getFileAddress(directory, std::string("newfx.dat"));
	std::ofstream force_x_file;
	force_x_file.open(force_x_file_address);

	std::string force_y_file_address = getFileAddress(directory, std::string("newfy.dat"));
	std::ofstream force_y_file;
	force_x_file.open(force_y_file_address);

	std::string force_z_file_address = getFileAddress(directory, std::string("newfz.dat"));
	std::ofstream force_z_file;
	force_z_file.open(force_z_file_address);

//Initialization.......................................
	double E_min, delE, E_max, T;
	int equil_steps, total_steps, m, TOTAL_ATOMS;

	input_file >> read_string >> TOTAL_ATOMS;
	input_file >> read_string >> equil_steps;
	input_file >> read_string >> total_steps;
	input_file >> read_string >> T;
	input_file >> read_string >> m;
	input_file >> read_string >> E_min;
	input_file >> read_string >> delE;
	input_file >> read_string >> E_max;

	int bins = std::ceil(((E_max-E_min)/delE) - 1) + 1;
	Vec_d E_bins; E_bins.resize(bins); E_bins.setZero();
	Vec_d E_tot_hist; E_tot_hist.resize(bins); E_tot_hist.setZero();
	Vec_d E_Prob_ref; E_Prob_ref.resize(bins); E_Prob_ref.setZero();

	for(int i = 0; i < bins; i++){
		E_bins[i] = E_min+i*delE;
		E_Prob_ref[i] = 1;
	}

//Equilibration phase..................................
	for (int i = 0; i < equil_steps; i++){
		double E;
		lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil"));
		lammps_file(LAMMPS_HANDLE, lammps_input_file_address.c_str());

		hist_file >> E;
		
		for(int j = 0; j < bins; j++){	
			if((E >= E_bins[j])&&(E < E_bins[j])){
			E_tot_hist[j]++;
			break;
			}
		}
	}

//CFES code............................................
	//Initializing...
	double old_E_factor = 0;
	Vec_d old_x; old_x.resize(TOTAL_ATOMS); old_x.setZero();
	Vec_d old_y; old_y.resize(TOTAL_ATOMS); old_y.setZero();
	Vec_d old_z; old_z.resize(TOTAL_ATOMS); old_z.setZero();

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes"));
	lammps_file(LAMMPS_HANDLE, lammps_input_file_address.c_str());

	hist_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	hist_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	xcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	xcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	ycoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	ycoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	zcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
	zcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"

	hist_file >> old_E_factor;
	for(int i = 0; i < TOTAL_ATOMS; i++){
		xcoord_file >> old_x[i];
		ycoord_file >> old_y[i];
		zcoord_file >> old_z[i];
	}
	xcoord_file.close();
	ycoord_file.close();
	zcoord_file.close();
	
	//Main loop...
	for (int i = 0; i < total_steps; i++){
		double E;
		Vec_d x; x.resize(TOTAL_ATOMS); x.setZero();
		Vec_d y; y.resize(TOTAL_ATOMS); y.setZero();
		Vec_d z; z.resize(TOTAL_ATOMS); z.setZero();

		hist_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		hist_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		xcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		xcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		ycoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		ycoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		zcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"
		zcoord_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"

		hist_file >> E;
		int j = 0;
		for(j = 0; j < bins; j++){	
			if((E >= E_bins[j])&&(E < E_bins[j])){
       			E_tot_hist[j]++;
			break;
			}
		}
		Vec_d E_Prob = E_tot_hist/E_tot_hist.sum(); 
	 	double E_factor = -(k_B*T)*std::log(1 - std::pow((E_Prob[j]/E_Prob_ref[j]), (m-1)));
      	        double dE_factor = (E_factor - old_E_factor);

		for(int j = 0; j < TOTAL_ATOMS; j++){
			xcoord_file >> x[j];
			ycoord_file >> y[j];
			zcoord_file >> z[j];
			force_x_file << j+1 << "  " << -dE_factor/(TOTAL_ATOMS*(x[j]-old_x[j])) << std::endl;
			force_y_file << j+1 << "  " << -dE_factor/(TOTAL_ATOMS*(y[j]-old_y[j])) << std::endl;
			force_z_file << j+1 << "  " << -dE_factor/(TOTAL_ATOMS*(z[j]-old_z[j])) << std::endl;
		}

		lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes_addforce"));
		lammps_file(LAMMPS_HANDLE, lammps_input_file_address.c_str());

		old_E_factor = E_factor;
		old_x = x;
		old_y = y;
		old_z = z;
		hist_file.close();
		xcoord_file.close();
		ycoord_file.close();
		zcoord_file.close();

		lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes"));
		lammps_file(LAMMPS_HANDLE, lammps_input_file_address.c_str());
		}

		lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_unfix"));
		lammps_file(LAMMPS_HANDLE, lammps_input_file_address.c_str());

//......................................
		lammps_close(LAMMPS_HANDLE);
		MPI_Finalize();
}
