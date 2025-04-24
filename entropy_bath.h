#pragma once
#include "declarations.h"	
#include "cfes.h"

double get_C_Q(){

	if (m < 0){	
	std::cout << " Negative m is not allowed ! " << std::endl;
	throw std::exception();
	}

	Vec_d C_Q_array; C_Q_array.resize(bins); C_Q_array.setZero();
	
	if (m > 1){
	#pragma omp parallel for shared(E_Prob, E_Prob_ref, C_Q_array)
        for(int j = 0; j < bins; j++){
                C_Q_array[j] = E_Prob[j]*std::pow((E_Prob[j]/E_Prob_ref[j]), (m-1));
	}
	
	return C_Q_array.sum();
	}
	else if (m == 1)
	return 1.0;
	else 
	return 0.0;	
}


void reset_Delta_s(double& C_Q, double& zeta_S, double& zeta_S_shift){

	Delta_s = 1.0/(zeta_S-zeta_S_shift) + C_Q;
}

double get_Delta_s(){

	return Delta_s;
}

void entropy_bath_initialize(double& zeta_S, double& zeta_S_shift){

	//==================================================
	//Bath initialization
	//==================================================
	double C_Q = get_C_Q();
	reset_Delta_s(C_Q, zeta_S, zeta_S_shift);
}

double entropy_bath(int& step, double& zeta_S, double& zeta_S_shift, int& coupling_time){

	//==================================================
	//Calculate CFES constant
	//==================================================
	double C_Q = get_C_Q();

	if(((step+restart_from_step) > 0) && ((step+restart_from_step) % coupling_time == 0))
	reset_Delta_s(C_Q, zeta_S, zeta_S_shift);

	double new_zeta_S = 1.0/(Delta_s - C_Q) + zeta_S_shift;

	if((m > 0) && (debug) && (mpi_id == 0) && (new_zeta_S < 0))
	std::cout << " zeta_S " << new_zeta_S << " is negative at step " << step << " ! " << std::endl;	

	return new_zeta_S;
}
