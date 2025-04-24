#pragma once
#include "declarations.h"	
#include "cfes.h"

void get_inputs(){

	//==================================================
	//Read inputs
	//==================================================

	input_file >> read_string >> read_string;
	if(read_string == std::string("yes"))
	debug = true;
	else debug = false;

	input_file >> read_string >> read_string;
	if(read_string == std::string("yes")){
	restart = true;
	input_file >> read_string >> probability_input_file_address;
	input_file >> read_string >> restart_from_step;
	}
	else {
	restart = false;
	restart_from_step = 0;
	}

	input_file >> read_string >> atom_index_file_address;
	input_file >> read_string >> total_steps;
	input_file >> read_string >> equil_steps;
	input_file >> read_string >> zeta_equil_steps;
	input_file >> read_string >> output_steps;
	input_file >> read_string >> min_prob_weight;
	input_file >> read_string >> max_prob_weight;
	input_file >> read_string >> E_min;
	input_file >> read_string >> E_max;
	input_file >> read_string >> binwidth;
	input_file >> read_string >> tau_equil;
	input_file >> read_string >> tau;
	input_file >> read_string >> zeta_S_ref;
	input_file >> read_string >> zeta_S_ref_shift;
	input_file >> read_string >> max_dE;
	input_file >> read_string >> k_b;
	input_file >> read_string >> T;
	input_file >> read_string >> m;
	input_file >> read_string >> mode;

	//==================================================
	//Read atom indices
	//==================================================

	std::ifstream atom_index_file(atom_index_file_address.c_str());
	atom_index_file >> read_string >> TOTAL_ATOMS;

	Atom_ID.resize(TOTAL_ATOMS); Atom_ID.setZero();

	for(int i = 0; i < TOTAL_ATOMS; i++){
		double index;
		atom_index_file >> index;
		Atom_ID[i] = index;
	}

	//==================================================
	//Set energy bins
	//==================================================

	bins = std::floor((E_max - E_min)/binwidth) + 1;
	E_bins.resize(bins); E_bins.setZero();
	E_tot_hist_equil.resize(bins); E_tot_hist_equil.setZero();
	E_tot_hist.resize(bins); E_tot_hist.setZero();
	E_Prob.resize(bins); E_Prob.setZero();
	E_Prob_equil.resize(bins); E_Prob_equil.setZero();
	E_Prob_ref.resize(bins); E_Prob_ref.setZero();
	E_tot_prodution_hist.resize(bins); E_tot_prodution_hist.setZero();

	for(int i = 0; i < bins; i++)
		E_bins[i] = E_min + i*binwidth;

	//==================================================
	//Set mode
	//==================================================
	
	if(mode == std::string("gaussian")){
	input_file >> read_string >> E_mean_incr_percentage;
	input_file >> read_string >> E_sigma_incr_percentage;
	input_file >> read_string >> E_norm_incr_percentage;
	input_file >> read_string >> sigma_shift_parameter; 
	}
	if(mode == std::string("constant"))
	input_file >> read_string >> sigma_shift_parameter; 
	if(mode == std::string("file")){//Always use file mode for restart !
	input_file >> read_string >> ref_prob_input_file_address;
	std::ifstream ref_prob_input_file(ref_prob_input_file_address.c_str());
	double energy;

	for(int i = 0; i < bins; i++)
		ref_prob_input_file >> energy >> E_Prob_ref[i];
	}
	if( (restart) && (mode != std::string("file")) ){
	std::cout << " Always use 'file' mode for restart ! " << std::endl; 
	throw std::exception();
	}

	//==================================================
	//Read restart probabilities (normalized)
	//==================================================

	if(min_prob_weight > max_prob_weight){
	std::cout << " min_prob_weight > max_prob_weight ! " << std::endl;
	throw std::exception();
	}

	if(restart){
	std::ifstream probability_input_file(probability_input_file_address.c_str());
	double energy;

	for(int i = 0; i < bins ; i++){
		probability_input_file >> energy >> E_Prob[i];
		E_tot_hist[i] = std::floor(E_Prob[i] * min_prob_weight);
	}
	}

	//==================================================
	//Debug output
	//==================================================
	
	if(mpi_id == 0){
	std::cout << " debug 	    	"  << debug 		<< std::endl;
	if(restart)
	std::cout << " Restart step 	"  << restart_from_step << std::endl;
	std::cout << " TOTAL_ATOMS  	"  << TOTAL_ATOMS	<< std::endl;
	std::cout << " total_steps  	"  << total_steps   	<< std::endl;
	std::cout << " equil_steps  	"  << equil_steps   	<< std::endl;
	std::cout << " output_steps 	"  << output_steps  	<< std::endl;
	std::cout << " min_prob_weight 	"  << min_prob_weight  	<< std::endl;
	std::cout << " max_prob_weight 	"  << max_prob_weight  	<< std::endl;
	std::cout << " E_min  	    	"  << E_min  		<< std::endl;
	std::cout << " E_max	    	"  << E_max         	<< std::endl;
	std::cout << " binwidth     	"  << binwidth      	<< std::endl;
	std::cout << " bins         	"  << bins          	<< std::endl;
	std::cout << " tau_equil    	"  << tau_equil     	<< std::endl;
	std::cout << " tau	    	"  << tau	        << std::endl;
	std::cout << " zeta_S_ref   	"  << zeta_S_ref    	<< std::endl;
	std::cout << " zeta_S_ref_shift "  << zeta_S_ref_shift  << std::endl;
	std::cout << " max_dE	    	"  << max_dE		<< std::endl;
	std::cout << " k_b 	    	"  << k_b           	<< std::endl;
	std::cout << " T 	    	"  << T 	        << std::endl;
	std::cout << " m            	"  << m 	        << std::endl;
	std::cout << " mode         	"  << mode          	<< std::endl;
	}

	if((mode == std::string("gaussian")) && (mpi_id == 0)){
	std::cout << " E_mean_incr_percentage  "  << E_mean_incr_percentage  << std::endl;
	std::cout << " E_sigma_incr_percentage "  << E_sigma_incr_percentage << std::endl;
	std::cout << " E_norm_incr_percentage  "  << E_norm_incr_percentage  << std::endl;
	std::cout << " sigma_shift_parameter   "  << sigma_shift_parameter   << std::endl;
	}

	if((mode == std::string("constant")) && (mpi_id == 0))
	std::cout << " sigma_shift_parameter   "  << sigma_shift_parameter   << std::endl;
}
