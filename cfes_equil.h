#pragma once
#include "declarations.h"	
#include "cfes.h"

void equilibrate(){	

	//==================================================
	//Equilibration
	//==================================================

	if(mpi_id == 0)
	std::cout << "Starting normal MD equilibration ... " << std::endl;

	int new_equil_steps = std::ceil((double)equil_steps/(double)tau_equil);
	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");

	for(int i = 0; i < new_equil_steps; i++){
		if(mpi_id == 0)
		print_progress_bar<int>(i, new_equil_steps, 10);

		std::string command = std::string("run ")+std::to_string(tau_equil)+std::string(" pre no post no");
		lmp->input->one(command.c_str());
		
		double *PE = (double *) lammps_extract_variable(lmp,"TrgtE","TrgtGrp");
		double E = *PE;
		int E_id = presorted_search_array<double>(E, E_min, E_max, binwidth);
		if(!std::isnan(E_id))
		E_tot_hist_equil[E_id]++;
	}

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());

	E_Prob_equil = E_tot_hist_equil/E_tot_hist_equil.sum(); 
	get_gaussian(E_bins, E_Prob_equil, E_mean_equil, E_sigma_equil, E_norm_equil);
	E_Prob = E_Prob_equil;
	E_tot_hist = E_tot_hist_equil;

	for(int i = 0; i < bins; i++){
		E_tot_hist[i] = std::floor(E_Prob[i]*min_prob_weight);
		if(mpi_id == 0)
		probability_equil_file << std::setprecision(8) 
			               << E_bins[i] << "  " << E_Prob[i] << std::endl;
	}
	if(mpi_id == 0)
	probability_equil_file << "\n\n" << std::endl;
}

void zeta_increment(){
	//==================================================
	//CFES
	//==================================================

	if(mpi_id == 0)
	std::cout << " CFES equilibration: increasing zeta_S from 0 to " << zeta_S_ref << std::endl;

	E_factor = 0;
	dE_factor_orig = 0;
	dE_factor = 0;
	P_factor = 0;
	old_E_factor = 0;

	Mat_d old_xyz; old_xyz.resize(TOTAL_ATOMS, 3); old_xyz.setZero();

	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes_equil"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");

	//==================================================
	//Main loop
	//==================================================

	double new_zeta_S_ref = 0;
	double new_zeta_S_ref_shift = 0;
	entropy_bath_initialize(new_zeta_S_ref, new_zeta_S_ref_shift);

	for (int i = 0; i < zeta_equil_steps; i++){
		
		if((i > 0) && (i % tau_equil == 0)){
		new_zeta_S_ref = ((double)i+1)*zeta_S_ref/(double)zeta_equil_steps;
		new_zeta_S_ref_shift = ((double)i+1)*zeta_S_ref_shift/(double)zeta_equil_steps;
		entropy_bath_initialize(new_zeta_S_ref, new_zeta_S_ref_shift);
		}

		if(mpi_id == 0)
		print_progress_bar<int>(i, zeta_equil_steps, 10);

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

		int E_id = presorted_search_array<double>(E, E_min, E_max, binwidth);
		if(!std::isnan(E_id))
		E_tot_hist[E_id]++;
		E_Prob = E_tot_hist/E_tot_hist.sum(); 

		//==================================================
		//Apply CFES bias
		//==================================================

		double zeta_S = entropy_bath(i, new_zeta_S_ref, new_zeta_S_ref_shift, tau_equil);
		cfes_main(i, zeta_S, E_id);

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

				if(i > 0){
				if((dr_mod > 0) && (!std::isinf(dr_mod)) && (!std::isnan(dr_mod)))
				force[id][d] += (-dE_factor/(double(TOTAL_ATOMS)*dr_mod))*(dx/dr_mod);
				else{	
				if((debug) && (mpi_id == 0))
				std::cout << " WARNING ! particle did not move at step " << i << " ! " 
					  << " dr_mod = " << dr_mod << std::endl;
				throw std::runtime_error(std::string("CFES simulation failed!"));
				}
				}
				old_xyz(id, d) = xyz[id][d];
			}
		}

	}
	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes_equil.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());
	
	#pragma omp parallel for shared(E_bins, E_Prob, E_tot_hist)
	for(int i = 0; i < bins; i++){
		if(mpi_id == 0)
		probability_equil_file << std::setprecision(8) 
			               << E_bins[i] << "  " << E_Prob[i] << std::endl;
	}
	if(mpi_id == 0)
	probability_equil_file << "\n\n" << std::endl;
}
