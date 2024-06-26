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
	double k_b, T;
	int bins, output_steps, equil_steps, total_steps, tau, m, TOTAL_ATOMS;
	std::string ref_Prob_file_address, atom_index_file_address;

	input_file >> read_string >> atom_index_file_address;
	input_file >> read_string >> total_steps;
	input_file >> read_string >> equil_steps;
	input_file >> read_string >> output_steps;
	input_file >> read_string >> ref_Prob_file_address;
	input_file >> read_string >> tau;
	input_file >> read_string >> k_b;
	input_file >> read_string >> T;
	input_file >> read_string >> m;

	std::ifstream atom_index_file(atom_index_file_address.c_str());
	atom_index_file >> read_string >> TOTAL_ATOMS;

	Vec_i Atom_ID; Atom_ID.resize(TOTAL_ATOMS); Atom_ID.setZero();

	for(int i = 0; i < TOTAL_ATOMS; i++){
		double index;
		atom_index_file >> index;
		Atom_ID[i] = index;
	}

	std::ifstream ref_Prob_file(ref_Prob_file_address.c_str());
	ref_Prob_file >> read_string >> bins;

	Vec_d E_bins; E_bins.resize(bins); E_bins.setZero();
	Vec_d E_tot_hist; E_tot_hist.resize(bins); E_tot_hist.setZero();
	Vec_d E_Prob_ref; E_Prob_ref.resize(bins); E_Prob_ref.setZero();

	for(int i = 0; i < bins; i++){
		double e, p;
		ref_Prob_file >> e >> p;
		E_bins[i] = e;
		E_Prob_ref[i] = p;
	}

	std::cout << " TOTAL_ATOMS  "  << TOTAL_ATOMS  << std::endl;
	std::cout << " total_steps  "  << total_steps  << std::endl;
	std::cout << " equil_steps  "  << equil_steps  << std::endl;
	std::cout << " output_steps "  << output_steps << std::endl;
	std::cout << " tau 	    "  << tau 	       << std::endl;
	std::cout << " T 	    "  << T 	       << std::endl;
	std::cout << " m            "  << m 	       << std::endl;

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
		
		double *PE = (double *) lammps_extract_variable(lmp,"TrgtE","TrgtGrp");
		double E = *PE;
		
		#pragma omp parallel for shared(E, E_tot_hist, E_bins)
		for(int j = 0; j < (bins-1); j++){	
			if((E >= E_bins[j]) && (E < E_bins[j+1])){
       			E_tot_hist[j]++;
			}
		}
	}

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());

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
		double *PE = (double *) lammps_extract_variable(lmp,"TrgtE","TrgtGrp");
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
	
		if((E >= E_bins[0]) && (E <= E_bins[bins-1])){
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
			       << E_bins[0] << " - " << E_bins[bins-1] << std::endl;

		#pragma omp parallel for shared(xyz, old_xyz, force)
		for(int i = 0; i < TOTAL_ATOMS; i++){
			int id = Atom_ID[i]-1;

			for(int d = 0; d < 3; d++){
				double dx = xyz[id][d] - old_xyz(id, d);
				if((i > 0) && (dx > 0) && (i % tau == 0))
				force[id][d] += -dE_factor/(TOTAL_ATOMS*dx);
				old_xyz(id, d) = xyz[id][d];
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
		probability_file << "\n" << std::endl;
		}
		}

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());

//......................................
	delete lmp;
	MPI_Finalize();
}
