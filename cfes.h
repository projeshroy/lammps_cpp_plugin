#pragma once
#include "declarations.h"	

	//Variables and Constants.............................
	void   *LAMMPS_HANDLE;
	LAMMPS *lmp;
	const double k_b_CGS = 0.001987204259; // Units: kcal/mol.K
	const double k_b_eV = 0.000086173303; // Units: eV/K
	const double PI = 3.14159265359;

	std::ifstream input_file;
	std::ofstream e_factor_file;
	std::ofstream e_factor_bath_file;
	std::ofstream probability_file;
	std::ofstream probability_final_file;
	std::ofstream probability_equil_file;
	std::ofstream ref_prob_file;
	std::ifstream atom_index_file;
 
	double  k_b, T, E, E_factor, P_factor, dE_factor_orig, dE_factor, old_E_factor, E_min, E_max, binwidth, max_dE,
		E_mean_equil, E_sigma_equil, E_norm_equil,
		zeta_S_ref, zeta_S_ref_shift, Delta_s,
       		E_mean_incr_percentage, E_sigma_incr_percentage, E_norm_incr_percentage, 
		sigma_shift_parameter, min_prob_weight, max_prob_weight;

	int mpi_total, mpi_id, output_steps, equil_steps, zeta_equil_steps, total_steps, tau, m, 
	    TOTAL_ATOMS, tau_equil, bins, sigma_shift_parameter_steps, restart_from_step;

	bool debug, restart;

	std::string directory, probability_input_file_address, atom_index_file_address, 
		    lammps_restart_input_file_address, ref_prob_input_file_address, mode;

	Vec_i Atom_ID; 	

	Vec_d E_bins, E_tot_hist_equil, E_tot_hist, E_Prob, E_Prob_equil, E_Prob_ref, E_tot_prodution_hist;

