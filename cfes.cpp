#include "cfes.h"
#include "Start.h"
#include "get_ref_probs.h"
#include "entropy_bath.h"
#include "cfes_main.h"
#include "cfes_equil.h"
#include "get_inputs.h"

int main(int argc, char **argv)
{
	//==================================================
	//LAMMPS initialization
	//==================================================

	directory = std::string("./");
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_total);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
	Start(argc, argv);

	//==================================================
	//Files
	//==================================================

	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps"));
	lmp->input->file(lammps_input_file_address.c_str());
		
	std::string input_file_address = getFileAddress(directory, std::string("input.in"));
	input_file.open(input_file_address.c_str());

	if(mpi_id == 0){
	std::string e_factor_file_address = getFileAddress(directory, std::string("E_factor.dat"));
	e_factor_file.open(e_factor_file_address, std::ios::out | std::ios::app);

	std::string e_factor_bath_file_address = getFileAddress(directory, std::string("E_factor_bath.dat"));
	e_factor_bath_file.open(e_factor_bath_file_address, std::ios::out | std::ios::app);

	std::string probability_equil_file_address = getFileAddress(directory, std::string("Probability_equil.dat"));
	probability_equil_file.open(probability_equil_file_address, std::ios::out | std::ios::app);

	std::string probability_file_address = getFileAddress(directory, std::string("Probability.dat"));
	probability_file.open(probability_file_address, std::ios::out | std::ios::app);

	std::string probability_final_file_address = getFileAddress(directory, std::string("Probability_final.dat"));
	probability_final_file.open(probability_final_file_address, std::ios::out | std::ios::app);

	std::string ref_prob_file_address = getFileAddress(directory, std::string("ref_prob.dat"));
	ref_prob_file.open(ref_prob_file_address, std::ios::out | std::ios::app);
	}

	//==================================================
	//Initialization
	//==================================================

	int step = 0;
	get_inputs();

	if(!restart){
	equilibrate();
	get_ref_probs(step);
	zeta_increment();
	}

	entropy_bath_initialize(zeta_S_ref, zeta_S_ref_shift);
	if((debug) && (mpi_id == 0))
	std::cout << " Initial Delta_s " << get_Delta_s() << " zeta_S " << zeta_S_ref << std::endl;

	//==================================================
	//CFES
	//==================================================

	if(mpi_id == 0)
	std::cout << " Initializing CFES ... " << std::endl;

	E_factor = 0;
	dE_factor_orig = 0;
	dE_factor = 0;
	P_factor = 0;
	old_E_factor = 0;

	Mat_d old_xyz; old_xyz.resize(TOTAL_ATOMS, 3); old_xyz.setZero();

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes"));
	lmp->input->file(lammps_input_file_address.c_str());
	if(restart){
	std::string command = std::string("reset_timestep  ") + std::to_string(restart_from_step);
	lmp->input->one(command.c_str());
	}
	lmp->input->one("run 0");

	if((!restart) && (mpi_id == 0)){
	e_factor_file << "# 1.Step  2.Energy  3.P_factor  4.E_factor  5.E_factor/(k_b*T)  6.dE_factor_orig   7.dE_factor " << std::endl;
	e_factor_bath_file << "# 1.Step  2.Delta_s  3.C_Q  4.zeta_S " << std::endl;
	}

	//==================================================
	//Main loop
	//==================================================

	for (step = restart_from_step; step < total_steps; step++){
		if(mpi_id == 0)
		print_progress_bar<int>(step, total_steps, 10);

		std::string command = std::string("run 1 pre no post no");
		lmp->input->one(command.c_str());
		
		//==================================================
		//LAMMPS data
		//==================================================

		double **xyz = lmp->atom->x;
		double **force = lmp->atom->f;
		double *PE = (double *) lammps_extract_variable(lmp,"TrgtE","TrgtGrp");
		E = *PE;

		//==================================================
		//Update histogram and bath
		//==================================================
		if(E_tot_hist.sum() >= max_prob_weight){
		#pragma omp parallel for shared(E_bins, E_Prob, E_tot_hist)
		for(int j = 0; j < bins; j++)
			E_tot_hist[j] = std::floor(E_Prob[j]*min_prob_weight);
		}

		double E_per_atom = E/(double)TOTAL_ATOMS;
		int E_id = presorted_search_array<double>(E, E_min, E_max, binwidth);
		if(!std::isnan(E_id)){
		E_tot_hist[E_id]++;
		E_tot_prodution_hist[E_id]++;
		}
		E_Prob = E_tot_hist/E_tot_hist.sum(); 
		
		//==================================================
		//Apply CFES bias
		//==================================================

		double zeta_S = entropy_bath(step, zeta_S_ref, zeta_S_ref_shift, tau);
		cfes_main(step, zeta_S, E_id);
	
		//==================================================
		//Update force
		//==================================================

		#pragma omp parallel for shared(xyz, old_xyz, force, dE_factor)
		for(int a = 0; a < TOTAL_ATOMS; a++){
			int id = Atom_ID[a]-1;
			double dr_mod = 0;

			for(int d = 0; d < 3; d++)
				dr_mod += pow((xyz[id][d] - old_xyz(id, d)), 2);

			dr_mod = sqrt(dr_mod);

			for(int d = 0; d < 3; d++){
				double dx = xyz[id][d] - old_xyz(id, d);

				if(step > 0){
				if((dr_mod > 0) && (!std::isinf(dr_mod)) && (!std::isnan(dr_mod)))
				force[id][d] += (-dE_factor/(double(TOTAL_ATOMS)*dr_mod))*(dx/dr_mod);
				else{	
				if((debug) && (mpi_id == 0))
				std::cout << " WARNING ! particle did not move at step " << step << " ! " 
					  << " dr_mod = " << dr_mod << std::endl;
				throw std::runtime_error(std::string("CFES simulation failed!"));
				}
				}
				old_xyz(id, d) = xyz[id][d];
			}
		}	
		//==================================================
		//Output
		//==================================================

		if(mpi_id == 0){
		if((step+restart_from_step) % output_steps == 0){
		double C_Q = get_C_Q();
		e_factor_file << std::setprecision(8) 
			      << step+restart_from_step << "  " << E << "  " << P_factor << "  " 
			      << E_factor << "  " << E_factor/(k_b*T) << "  " 
			      << dE_factor_orig << "  " << dE_factor << std::endl;

		e_factor_bath_file << std::setprecision(8) 
				   << step+restart_from_step << "  " << Delta_s << "  " 
				   << C_Q << "  " << zeta_S << std::endl;

		for(int j = 0; j < bins; j++)
			probability_file << std::setprecision(8) 
				         << E_bins[j] << "  " << E_Prob[j] << std::endl;
			
		probability_file << "\n" << std::endl;
		}}
		//==================================================
		}
	
	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());

	for(int j = 0; j < bins; j++)
		probability_final_file  << std::setprecision(8) 
			         	<< E_bins[j] << "  " << E_tot_prodution_hist[j]/E_tot_prodution_hist.sum() << std::endl;
			
	probability_final_file << "\n" << std::endl;

//==================================================

	delete lmp;
	MPI_Finalize();
}
