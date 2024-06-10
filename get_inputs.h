#pragma once
#include "declarations.h"	
#include "cfes.h"

void get_inputs(){

	input_file >> read_string >> atom_index_file_address;
	input_file >> read_string >> total_steps;
	input_file >> read_string >> equil_steps;
	input_file >> read_string >> output_steps;
	input_file >> read_string >> E_min;
	input_file >> read_string >> E_max;
	input_file >> read_string >> binwidth;
	input_file >> read_string >> tau_equil;
	input_file >> read_string >> tau;
	input_file >> read_string >> C_Q_factor;
	input_file >> read_string >> max_dE;
	input_file >> read_string >> k_b;
	input_file >> read_string >> T;
	input_file >> read_string >> m;
	input_file >> read_string >> mode;
	
	if(mode == std::string("gaussian")){
	input_file >> read_string >> E_mean_incr_percentage;
	input_file >> read_string >> E_sigma_incr_percentage;
	input_file >> read_string >> E_norm_incr_percentage;
	input_file >> read_string >> sigma_shift_parameter; 
	}
	if(mode == std::string("constant"))
	input_file >> read_string >> sigma_shift_parameter; 
	if(mode == std::string("dynamic")){
	input_file >> read_string >> sigma_shift_parameter_target;
	input_file >> read_string >> sigma_shift_parameter_increment;
 	input_file >> read_string >> sigma_shift_parameter_steps;
	}

	std::ifstream atom_index_file(atom_index_file_address.c_str());
	atom_index_file >> read_string >> TOTAL_ATOMS;

	Atom_ID.resize(TOTAL_ATOMS); Atom_ID.setZero();

	for(int i = 0; i < TOTAL_ATOMS; i++){
		double index;
		atom_index_file >> index;
		Atom_ID[i] = index;
	}

	bins = std::floor((E_max - E_min)/binwidth) + 1;
	E_bins.resize(bins); E_bins.setZero();
	E_tot_hist_equil.resize(bins); E_tot_hist_equil.setZero();
	E_tot_hist.resize(bins); E_tot_hist.setZero();
	E_Prob.resize(bins); E_Prob.setZero();
	E_Prob_equil.resize(bins); E_Prob_equil.setZero();
	E_Prob_ref.resize(bins); E_Prob_ref.setZero();

	for(int i = 0; i < bins; i++)
		E_bins[i] = E_min + i*binwidth;

	std::cout << " TOTAL_ATOMS  "  << TOTAL_ATOMS	<< std::endl;
	std::cout << " total_steps  "  << total_steps   << std::endl;
	std::cout << " equil_steps  "  << equil_steps   << std::endl;
	std::cout << " output_steps "  << output_steps  << std::endl;
	std::cout << " E_min  	    "  << E_min  	<< std::endl;
	std::cout << " E_max	    "  << E_max         << std::endl;
	std::cout << " binwidth     "  << binwidth      << std::endl;
	std::cout << " bins         "  << bins          << std::endl;
	std::cout << " tau_equil    "  << tau_equil     << std::endl;
	std::cout << " tau	    "  << tau	        << std::endl;
	std::cout << " C_Q_factor   "  << C_Q_factor    << std::endl;
	std::cout << " max_dE	    "  << max_dE	<< std::endl;
	std::cout << " k_b 	    "  << k_b           << std::endl;
	std::cout << " T 	    "  << T 	        << std::endl;
	std::cout << " m            "  << m 	        << std::endl;
	std::cout << " mode         "  << mode          << std::endl;

	//==================================================
	//Mode == Gaussian
	//==================================================

	if(mode == std::string("gaussian")){
	std::cout << " E_mean_incr_percentage  "  << E_mean_incr_percentage  << std::endl;
	std::cout << " E_sigma_incr_percentage "  << E_sigma_incr_percentage << std::endl;
	std::cout << " E_norm_incr_percentage  "  << E_norm_incr_percentage  << std::endl;
	std::cout << " sigma_shift_parameter   "  << sigma_shift_parameter   << std::endl;
	}

	//==================================================
	//Mode == Constant
	//==================================================

	if(mode == std::string("constant"))
	std::cout << " sigma_shift_parameter   "  << sigma_shift_parameter   << std::endl;

	//==================================================
	//Mode == Dynamic_constant
	//==================================================

	if(mode == std::string("dynamic")){
	std::cout << " sigma_shift_parameter_target     "  << sigma_shift_parameter_target      << std::endl;
	std::cout << " sigma_shift_parameter_increment  "  << sigma_shift_parameter_increment   << std::endl;
	std::cout << " sigma_shift_parameter_steps 	"  << sigma_shift_parameter_steps   	<< std::endl;
	}
}
