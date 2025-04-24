#include "Start.h"

int main(int argc, char **argv)
{
//LAMMPS initialization................................
	std::string directory = argv[1];
	Start(argc, argv);

//Files.................................................
	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps"));
	lmp->input->file(lammps_input_file_address.c_str());

	std::string input_file_address = getFileAddress(directory, std::string("input.in"));
	std::ifstream input_file(input_file_address.c_str());

	std::string e_factor_file_address = getFileAddress(directory, std::string("E_factor.dat"));
	std::ofstream e_factor_file;
	e_factor_file.open(e_factor_file_address);

	std::string force_x_file_address = getFileAddress(directory, std::string("newfx.dat"));
	std::ofstream force_x_file;
	force_x_file.open(force_x_file_address);

	std::string force_y_file_address = getFileAddress(directory, std::string("newfy.dat"));
	std::ofstream force_y_file;
	force_y_file.open(force_y_file_address);

	std::string force_z_file_address = getFileAddress(directory, std::string("newfz.dat"));
	std::ofstream force_z_file;
	force_z_file.open(force_z_file_address);

//Initialization.......................................
	double E_min, delE, E_max, T;
	int cleanup_step, equil_steps, total_steps, m, TOTAL_ATOMS;

	input_file >> read_string >> TOTAL_ATOMS;
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

	std::cout << " TOTAL_ATOMS "  << TOTAL_ATOMS  << std::endl;
	std::cout << " total_steps "  << total_steps  << std::endl;
	std::cout << " T 	   "  << T 	      << std::endl;
	std::cout << " m           "  << m 	      << std::endl;
	std::cout << " E_min       "  << E_min 	      << std::endl;
	std::cout << " delE        "  << delE 	      << std::endl;
	std::cout << " E_max       "  << E_max        << std::endl;

//Equilibration phase..................................
	std::cout << "Starting Equilibration ... " << std::endl;

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil"));
	lmp->input->file(lammps_input_file_address.c_str());

	std::string hist_file_address = getFileAddress(directory, std::string("hist.dat"));
	std::ifstream hist_equil_file(hist_file_address.c_str());
	hist_equil_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//ignore "#"

	while (true) {
	double E;
	hist_equil_file >> E;
	
	#pragma omp parallel for shared(E, E_bins, E_tot_hist)
	for(int j = 0; j < (bins-1); j++){	
		if((E >= E_bins[j])&&(E < E_bins[j+1])){
     		E_tot_hist[j]++;
		}
	}
	if(hist_equil_file.eof()) break;
	}

	hist_equil_file.close();

//CFES code............................................	
	std::cout << " Initializing CFES ... " << std::endl;

	double old_E_factor = 0;
	double* old_x;
	double* old_y;
	double* old_z;

	for(int i = 0; i < TOTAL_ATOMS; i++){
		force_x_file << i+1 << "  " << 0.0 << std::endl;
		force_y_file << i+1 << "  " << 0.0 << std::endl;
		force_z_file << i+1 << "  " << 0.0 << std::endl;
	}

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");

	//Main loop...
	for (int i = 0; i < total_steps; i++){
		print_progress_bar<int>(i, total_steps, 10);
		lmp->input->one("run 1 pre no post no");
	
		//get coordinates
		double *x = (double *) lammps_extract_variable(lmp,"x","peptide");
		double *y = (double *) lammps_extract_variable(lmp,"y","peptide");
		double *z = (double *) lammps_extract_variable(lmp,"z","peptide");

//		auto xyz = lmp->atom->x;
		//get energy
		double *PE = (double *) lammps_extract_variable(lmp,"pepE","peptide");
		double E = *PE;
		double E_factor = 0;
		double dE_factor = 0;
		int E_id = 0;

		#pragma omp parallel for shared(E, E_tot_hist, E_bins, E_id)
		for(int j = 0; j < (bins-1); j++){	
			if((E >= E_bins[j])&&(E < E_bins[j+1])){
       			E_tot_hist[j]++;
			E_id = j;
			}
		}

		Vec_d E_Prob = E_tot_hist/E_tot_hist.sum(); 
	 	E_factor = -(k_B*T)*std::log(1 - std::pow((E_Prob[E_id]/E_Prob_ref[E_id]), (m-1)));	
		dE_factor = (E_factor - old_E_factor);
		e_factor_file << i << "  " << E << "  " << E_factor << "  " << dE_factor << std::endl;

		for(int j = 0; j < TOTAL_ATOMS; j++){
			if(i > 0){
                        force_x_file << j+1 << "  " << -dE_factor/(TOTAL_ATOMS*(x[j]-old_x[j])) << std::endl;
 		        force_y_file << j+1 << "  " << -dE_factor/(TOTAL_ATOMS*(y[j]-old_y[j])) << std::endl;
		        force_z_file << j+1 << "  " << -dE_factor/(TOTAL_ATOMS*(z[j]-old_z[j])) << std::endl;
		}
		}

		old_E_factor = E_factor;
		old_x = x;
		old_y = y;
		old_z = z;
		}

//......................................
	delete lmp;
	MPI_Finalize();
}
