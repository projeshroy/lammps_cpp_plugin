#pragma once
#include "declarations.h"	
#include "cfes.h"

void equilibrate(){	

	//==================================================
	//Equilibration
	//==================================================

	std::cout << "Starting Equilibration ... " << std::endl;

	int new_equil_steps = std::ceil((double)equil_steps/(double)tau_equil);
	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");

	for(int i = 0; i < new_equil_steps; i++){
		print_progress_bar<int>(i, new_equil_steps, 10);
		std::string command = std::string("run ")+std::to_string(tau_equil)+std::string(" pre no post no");
		lmp->input->one(command.c_str());
		
		double *PE = (double *) lammps_extract_variable(lmp,"TrgtE","TrgtGrp");
		double E = *PE;
		
		#pragma omp parallel for shared(E, E_tot_hist, E_bins)
		for(int j = 0; j < (bins-1); j++){	
			if((E >= E_bins[j]) && (E < E_bins[j+1])){
       			E_tot_hist_equil[j]++;
			}
		}
	}

	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_equil.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());

	E_Prob_equil = E_tot_hist_equil/E_tot_hist_equil.sum(); 
	get_gaussian(E_bins, E_Prob_equil, E_mean_equil, E_sigma_equil, E_norm_equil);
	E_Prob = E_Prob_equil;
	for(int i = 0; i < bins; i++)
		E_tot_hist[i] = std::floor(E_Prob[i]*tau_equil);

}
