#include "Start.h"

int main(int argc, char **argv)
{
//LAMMPS initialization................................
	std::string directory = std::string("./");
	Start(argc, argv);

//Files.................................................
	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps"));
	lmp->input->file(lammps_input_file_address.c_str());

	std::string input_file_address = getFileAddress(directory, std::string("input.in"));
	std::ifstream input_file(input_file_address.c_str());

	std::string e_factor_file_address = getFileAddress(directory, std::string("E_factor.dat"));
	std::ofstream e_factor_file;
	e_factor_file.open(e_factor_file_address);

	std::string probability_file_address = getFileAddress(directory, std::string("Probability.dat"));
	std::ofstream probability_file;
	probability_file.open(probability_file_address);

//Initialization.......................................
	double E_min, delE, E_max, T;
	int output_steps, equil_steps, total_steps, tau, m, TOTAL_ATOMS;

	input_file >> read_string >> TOTAL_ATOMS;
	input_file >> read_string >> total_steps;
	input_file >> read_string >> equil_steps;
	input_file >> read_string >> output_steps;
	input_file >> read_string >> tau;
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

	std::cout << " TOTAL_ATOMS  "  << TOTAL_ATOMS  << std::endl;
	std::cout << " total_steps  "  << total_steps  << std::endl;
	std::cout << " equil_steps  "  << equil_steps  << std::endl;
	std::cout << " output_steps "  << output_steps << std::endl;
	std::cout << " tau 	    "  << tau 	       << std::endl;
	std::cout << " T 	    "  << T 	       << std::endl;
	std::cout << " m            "  << m 	       << std::endl;
	std::cout << " E_min        "  << E_min        << std::endl;
	std::cout << " delE         "  << delE 	       << std::endl;
	std::cout << " E_max        "  << E_max        << std::endl;

	int new_equil_steps = std::ceil((double)equil_steps/(double)tau);

//Equilibration phase..................................
	std::cout << "Starting Equilibration ... " << std::endl;

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");

	for(int i = 0; i < new_equil_steps; i++){
		print_progress_bar<int>(i, new_equil_steps, 10);
		std::string command = std::string("run ")+std::to_string(tau)+std::string(" pre no post no");
		lmp->input->one(command.c_str());
		
		double *PE = (double *) lammps_extract_variable(lmp,"pepE","peptide");
		double E = *PE;

		#pragma omp parallel for shared(E, E_tot_hist, E_bins)
		for(int j = 0; j < (bins-1); j++){	
			if((E >= E_bins[j]) && (E < E_bins[j+1])){
       			E_tot_hist[j]++;
			}
		}
	}

	lmp->input->one("unfix   pepmd");
	lmp->input->one("unfix   mom");
	lmp->input->one("unfix   file");
	lmp->input->one("unfix   constraint");
	lmp->input->one("undump  traj");

//CFES code............................................	
	std::cout << " Initializing CFES ... " << std::endl;

	double old_E_factor = 0;
	Mat_d old_xyz; old_xyz.resize(TOTAL_ATOMS, 3); old_xyz.setZero();
	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");
	e_factor_file << "# 1.Step  2.Energy  3.Probability  4.E_factor  5.E_factor/(k_B*T)  6.dE_factor " << std::endl;

	//Main loop...
	for (int i = 0; i < total_steps; i++){
		print_progress_bar<int>(i, total_steps, 10);
		std::string command = std::string("run 1 pre no post no");
		lmp->input->one(command.c_str());
		
		//..................................................	
		//LAMMPS data
		//..................................................
		double **xyz = lmp->atom->x;
		double **force = lmp->atom->f;
		double *PE = (double *) lammps_extract_variable(lmp,"pepE","peptide");
		double E = *PE;

		//..................................................	
		//Update histogram
		//..................................................
		int E_id = 0;

		if((i > 0) && (i % tau == 0)){
		#pragma omp parallel for shared(E, E_tot_hist, E_bins, E_id)
		for(int j = 0; j < (bins-1); j++){	
			if((E >= E_bins[j])&&(E < E_bins[j+1])){
       			E_tot_hist[j]++;
			E_id = j;
			}
		}
		}
		Vec_d E_Prob = E_tot_hist/E_tot_hist.sum(); 

		//..................................................	
		//CFES bias
		//..................................................
		double E_factor = 0;
		double dE_factor = 0;
		double P_factor = 0;
	
		if((E >= E_min) && (E <= E_max)){
		P_factor = std::pow((E_Prob[E_id]/E_Prob_ref[E_id]), (m-1));
		if(P_factor < 1){
		E_factor = -(k_B*T)*std::log(1 - P_factor);
		dE_factor = (E_factor - old_E_factor);
		}
		else if(i % tau == 0) 
		     std::cout << " WARNING ! Probability " << E_Prob[E_id] 
		   	       << " at step " << i*tau 
			       << " is larger than the reference probability " << E_Prob_ref[E_id] << std::endl;	
		}
		else if(i % tau == 0)
		     std::cout << " WARNING ! Energy " << E 
			       << " is outside the range " 
			       << E_min << " - " << E_max << std::endl;

		#pragma omp parallel for shared(xyz, old_xyz, force)
		for(int i = 0; i < TOTAL_ATOMS; i++){
			for(int d = 0; d < 3; d++){
				double dx = xyz[i][d] - old_xyz(i, d);
				if((i > 0) && (dx > 0) && (i % tau == 0))
				force[i][d] += -dE_factor/(TOTAL_ATOMS*dx);
				old_xyz(i, d) = xyz[i][d];
			}
		}
		old_E_factor = E_factor;

		//..................................................	
		//Output
		//..................................................
		if(i % output_steps == 0){
		e_factor_file << i << "  " << E << "  " << E_Prob[E_id] << "  " << E_factor 
			      << "  " << E_factor/(k_B*T) << "  " << dE_factor << std::endl;
		for(int j = 0; j < bins; j++)
			probability_file << E_bins[j] << "  " << E_Prob[j] << std::endl;
		probability_file << "\n \n" << std::endl;
		}
		}

//......................................
	delete lmp;
	MPI_Finalize();
}
