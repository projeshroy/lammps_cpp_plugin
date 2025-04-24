#pragma once
#include "declarations.h"	
#include "cfes.h"

void set_C_Q_ref(){

	C_Q_ref = (1.0/(pow((C_Q - 1.0), 2))) + C_Q + C_Q_factor;
}

void entropy_bath_initialize(){
	//==================================================
	//Bath initialization
	//==================================================

	E_Prob_bath = E_Prob_equil;
	E_tot_hist_bath = E_tot_hist_equil;
	E_bath = E_mean_equil;
	P_factor_bath = 0;
	E_factor_bath = 0;
	lambda_S = 0;
	C_Q = 0;
	C_Q_ref = 0;

	if(m > 0){
	for(int j = 0; j < bins; j++)
		C_Q += E_Prob[j]*std::pow((E_Prob[j]/E_Prob_ref[j]), (m-1));

	set_C_Q_ref();
	if((debug) && (mpi_id == 0))
	std::cout << " Initial C_Q_ref " << C_Q_ref << std::endl;

	for(int j = 0; j < bins; j++)
		lambda_S += E_Prob_bath[j]/(C_Q_ref - C_Q);	
	}
}

void entropy_bath(int& step){

	//==================================================
	//Calculate CFES constant
	//==================================================

	C_Q = 0;
	for(int j = 0; j < bins; j++)
		C_Q += E_Prob[j]*std::pow((E_Prob[j]/E_Prob_ref[j]), (m-1));

	set_C_Q_ref();
	P_factor_bath = C_Q/C_Q_ref;
	E_factor_bath = -(TOTAL_ATOMS*k_b*T)*std::log(1.0 - P_factor_bath);
	E_bath = E_factor_bath + E_mean_equil; 

	//TOTAL_ATOMS and E_mean_equil are added here to scale it to E_bins only! 
	//It has no effect on the probability distribution.

	//==================================================
	//Update bath
	//==================================================

	int E_id = 0;
	#pragma omp parallel for shared(E_bath, E_tot_hist_bath, E_bins, E_id)
	for(int j = 0; j < (bins-1); j++){	
		if((E_bath >= E_bins[j]) && (E_bath < E_bins[j+1])){
		E_tot_hist_bath[j]++;
		E_id = j;
		}
	}
	E_Prob_bath = E_tot_hist_bath/E_tot_hist_bath.sum(); 

	lambda_S = 0;
	for(int j = 0; j < bins; j++)
		lambda_S += E_Prob_bath[j]/(C_Q_ref - C_Q);
	
}
