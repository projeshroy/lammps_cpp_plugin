bool cfes_main(int& step, bool& cfes_equil, double& zeta_S, int& E_id){
	//==================================================
	//CFES bias
	//==================================================
	
	E_factor = 0;
	dE_factor_orig = 0;
	dE_factor = 0;
	P_factor = 0;
	bool skip = false;

	if((E >= E_min) && (E <= E_max)){
	P_factor = (double)m*zeta_S*std::pow((E_Prob[E_id]/E_Prob_ref[E_id]), (m-1));
	//make sure zeta_S is a positive number!
	
	//==================================================
	//Actual run
	//==================================================

	if(P_factor < 1){
	E_factor = -((double)TOTAL_ATOMS*k_b*T)*std::log(1.0 - P_factor);
	dE_factor_orig = (E_factor - old_E_factor);
	dE_factor = dE_factor_orig;

	//Fluctuation prevention step.......................
	if(sqrt(dE_factor_orig*dE_factor_orig) > (max_dE*TOTAL_ATOMS*k_b*T)){
	if(cfes_equil){
	dE_factor = (dE_factor_orig/sqrt(dE_factor_orig*dE_factor_orig))*max_dE*TOTAL_ATOMS*k_b*T;	
	if((m > 0) && (debug) && (mpi_id == 0))
	std::cout << " Max dE reached at step " << step << " with zeta_S " << zeta_S << " ! replacing dE with " << dE_factor << std::endl;
	}
	else {
	skip = true;
	if((m > 0) && (debug) && (mpi_id == 0))
        std::cout << " Max dE reached at step " << step << " with zeta_S " << zeta_S << " ! No force added !" << std::endl;
	}
	}
	}

	//==================================================
	//Equilibration related steps
	//==================================================

	else {
	if((m > 0) && (debug) && (mpi_id == 0))
	std::cout << " WARNING ! Probability factor " << P_factor << " greater than 1 ! " << " at step " << step << " with zeta_S " << zeta_S 
		  << " Actual pobability: " << E_Prob[E_id] << " Reference probability: " << E_Prob_ref[E_id] << std::endl;
	
	if(cfes_equil){
	E_factor = ((double)TOTAL_ATOMS*k_b*T)*std::log(1.0 + P_factor);
	dE_factor_orig = (E_factor - old_E_factor);
	dE_factor = dE_factor_orig;

	if(sqrt(dE_factor_orig*dE_factor_orig) > (max_dE*TOTAL_ATOMS*k_b*T)){
	dE_factor = (dE_factor_orig/sqrt(dE_factor_orig*dE_factor_orig))*max_dE*TOTAL_ATOMS*k_b*T;	
	if((m > 0) && (debug) && (mpi_id == 0))
	std::cout << " Max dE reached at step " << step << " with zeta_S " << zeta_S << " ! replacing dE with " << dE_factor << std::endl;
	}
	if((m > 0) && (debug) && (mpi_id == 0))
	std::cout << " Zeta_S sign changed ! " << std::endl;
	}
	else 
	{
	skip = true;
	if((m > 0) && (debug) && (mpi_id == 0))
	std::cout << " No force added ! " << std::endl;	
	}
	}
	}
	//==================================================
	//Out of energy range
	//==================================================

	else {
	skip = true;
	if((m > 0) && (debug) && (mpi_id == 0))
	std::cout << " WARNING ! Energy " << E 
	          << " is outside the range " 
	          << E_bins[0] << " - " << E_bins[bins-1] << " at step " << step << " with zeta_S " << zeta_S << std::endl;
	}

	old_E_factor = E_factor;
	
	return skip;
}
