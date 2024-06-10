#pragma once
#include "declarations.h"	

	//Variables and Constants.............................
	void   *LAMMPS_HANDLE;
	LAMMPS *lmp;
	const double k_b_CGS = 0.001987204259; // Units: kcal/mol.K
	const double PI = 3.14159265359;

	std::ifstream input_file;
	std::ofstream e_factor_file;
	std::ofstream e_factor_bath_file;
	std::ofstream probability_file;
	std::ofstream probability_bath_file;
	std::ofstream ref_prob_file;
	std::ifstream atom_index_file;

	double  k_b, T, E, E_factor, P_factor, dE_factor_orig, dE_factor, old_E_factor, E_min, E_max, binwidth, max_dE,
		E_mean_equil, E_sigma_equil, E_norm_equil,
		E_bath, E_factor_bath, P_factor_bath, C_Q, lambda_S, C_Q_ref, C_Q_factor,
       		E_mean_incr_percentage, E_sigma_incr_percentage, E_norm_incr_percentage, 
		sigma_shift_parameter, 
		sigma_shift_parameter_target, sigma_shift_parameter_increment; 

	int output_steps, equil_steps, total_steps, tau, m, TOTAL_ATOMS, tau_equil, bins, sigma_shift_parameter_steps;

	std::string directory, atom_index_file_address, mode;

	Vec_i Atom_ID; 	
	Vec_d E_bins, E_tot_hist_equil, E_tot_hist, E_tot_hist_bath, E_Prob, E_Prob_bath, E_Prob_equil, E_Prob_ref;
