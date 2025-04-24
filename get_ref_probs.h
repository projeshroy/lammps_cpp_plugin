#pragma once
#include "declarations.h"	
#include "cfes.h"

void get_gaussian(Vec_d& E, Vec_d& Prob, double& mean, double& sigma, double& norm){

	//==================================================
	//Gaussian Approximation
	//==================================================

	mean = 0;
	sigma = 0;
	norm = 0;	

	for(int i = 0; i < E.rows(); i++)
		mean += E[i]*Prob[i];

	for(int i = 0; i < E.rows(); i++)
		sigma += (pow((E[i] - mean), 2))*Prob[i];

	sigma = sqrt(sigma);

	for(int i = 0; i < E.rows(); i++)
		norm += exp(-pow(((E[i] - mean)/(sqrt(2.0)*sigma)), 2));		

	norm = 1.0/norm;

//	norm = Prob.maxCoeff();
}

void get_ref_probs(int& step){

	//==================================================
	//Mode == Gaussian
	//==================================================

	if(mode == std::string("gaussian")){
	double E_mean, E_sigma, E_norm, P_shift;

	P_shift = E_norm_equil*exp(-pow(((sigma_shift_parameter*E_sigma_equil)/(sqrt(2.0)*E_sigma_equil)), 2));
	E_mean = E_mean_equil*(100+E_mean_incr_percentage)/100;
	E_sigma = E_sigma_equil*(100+E_sigma_incr_percentage)/100;
	E_norm = E_norm_equil*(100+E_norm_incr_percentage)/100 - P_shift;

	if(mpi_id == 0){
	std::cout << " E_mean_equil  " << E_mean_equil  << " E_mean  " << E_mean  << std::endl;
	std::cout << " E_sigma_equil " << E_sigma_equil << " E_sigma " << E_sigma << std::endl;
	std::cout << " E_norm_equil  " << E_norm_equil  << " E_norm  " << E_norm + P_shift  << std::endl;
	std::cout << " P_shift " << P_shift << std::endl;
	}

	if(mpi_id == 0)
	ref_prob_file << "# E   P" << std::endl;
	
	for(int i = 0; i < bins; i++){
		E_Prob_ref[i] = (E_norm*exp(-pow(((E_bins[i] - E_mean)/(sqrt(2.0)*E_sigma)), 2)) + P_shift);
		if(mpi_id == 0)
		ref_prob_file << E_bins[i] << "   " << E_Prob_ref[i] << std::endl;
	}
	}

	//==================================================
	//Mode == Constant
	//==================================================

	if(mode == std::string("constant")){
	double P_shift = E_norm_equil*exp(-pow(((sigma_shift_parameter*E_sigma_equil)/(sqrt(2.0)*E_sigma_equil)), 2));
	if(mpi_id == 0)
	std::cout << " P_shift " << P_shift << std::endl;

	for(int i = 0; i < bins; i++){
		E_Prob_ref[i] = P_shift;
		if(mpi_id == 0)	
		ref_prob_file << E_bins[i] << "   " << E_Prob_ref[i] << std::endl;
	}
	}
}

