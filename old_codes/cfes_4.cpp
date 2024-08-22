#include "cfes.h"
#include "Start.h"
#include "get_ref_probs.h"
#include "cfes_equil.h"
#include "get_inputs.h"
#include "entropy_bath.h"

int main(int argc, char **argv)
{
	//==================================================
	//LAMMPS initialization
	//==================================================

	directory = std::string("./");
	Start(argc, argv);

	//==================================================
	//Files
	//==================================================

	std::string lammps_input_file_address = getFileAddress(directory, std::string("in.lammps"));
	lmp->input->file(lammps_input_file_address.c_str());
		
	std::string input_file_address = getFileAddress(directory, std::string("input.in"));
	input_file.open(input_file_address.c_str());

	if(mpi_id == 0){
	std::string e_factor_file_address = getFileAddress(directory, std::string("E_factor.dat"));
	e_factor_file.open(e_factor_file_address);

	std::string e_factor_bath_file_address = getFileAddress(directory, std::string("E_factor_bath.dat"));
	e_factor_bath_file.open(e_factor_bath_file_address);

	std::string probability_file_address = getFileAddress(directory, std::string("Probability.dat"));
	probability_file.open(probability_file_address);

//	std::string probability_bath_file_address = getFileAddress(directory, std::string("Probability_bath.dat"));
//	probability_bath_file.open(probability_bath_file_address);

	std::string ref_prob_file_address = getFileAddress(directory, std::string("ref_prob.dat"));
	ref_prob_file.open(ref_prob_file_address);
	}

	//==================================================
	//Initialization
	//==================================================

	int step = 0;
	get_inputs();
	equilibrate();
	get_ref_probs(step);
	entropy_bath_initialize();

	//==================================================
	//CFES
	//==================================================

	if(mpi_id == 0)
	std::cout << " Initializing CFES ... " << std::endl;

	E_factor = 0;
	E_factor_bath = 0;
	dE_factor_orig = 0;
	dE_factor = 0;
	P_factor = 0;
	P_factor_bath = 0;
	old_E_factor = 0;

	Mat_d old_xyz; old_xyz.resize(TOTAL_ATOMS, 3); old_xyz.setZero();
	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes"));
	lmp->input->file(lammps_input_file_address.c_str());
	lmp->input->one("run 0");

	if(mpi_id == 0){
	e_factor_file << "# 1.Step  2.Energy  3.P_factor  4.E_factor  5.E_factor/(k_b*T)  6.dE_factor_orig   7.dE_factor " << std::endl;
//	e_factor_bath_file << "# 1.Step  2.Energy  3.P_factor  4.E_factor  5.E_factor/(k_b*T)  6.C_Q_ref  7.C_Q  8.zeta_S " << std::endl;
	e_factor_bath_file << "# 1.Step  2.C_Q_ref  3.C_Q  4.zeta_S " << std::endl;
	}

	//==================================================
	//Main loop
	//==================================================

	for (step = 0; step < total_steps; step++){
		if(mpi_id == 0)
		print_progress_bar<int>(step, total_steps, 10);

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

		int E_id = 0;

		#pragma omp parallel for shared(E, E_tot_hist, E_bins, E_id)
		for(int j = 0; j < (bins-1); j++){	
			if((E >= E_bins[j])&&(E < E_bins[j+1])){
       			E_tot_hist[j]++;
			E_id = j;
			}
		}
		E_Prob = E_tot_hist/E_tot_hist.sum(); 
		if((mode == std::string("dynamic")) && 
		   (step % sigma_shift_parameter_steps == 0))
		get_ref_probs(step);		
		
		if(m > 0)
		entropy_bath(step);

		//==================================================
		//CFES bias
		//==================================================
		
		E_factor = 0;
		dE_factor_orig = 0;
		dE_factor = 0;
		P_factor = 0;


		if( m > 0) {
		if((E >= E_bins[0]) && (E <= E_bins[bins-1])){
		P_factor = (double)m*zeta_S*std::pow((E_Prob[E_id]/E_Prob_ref[E_id]), (m-1));
		//make sure zeta_S is a positive number!
		
		if(P_factor < 1){
		E_factor = -(TOTAL_ATOMS*k_b*T)*std::log(1.0 - P_factor);
		dE_factor_orig = (E_factor - old_E_factor);
		dE_factor = dE_factor_orig;

		if(sqrt(dE_factor_orig*dE_factor_orig) > (max_dE*TOTAL_ATOMS*k_b*T))
		dE_factor = (dE_factor_orig/sqrt(dE_factor_orig*dE_factor_orig))*max_dE*TOTAL_ATOMS*k_b*T;
		}
		else {
		E_factor = (TOTAL_ATOMS*k_b*T)*std::log(1.0 + P_factor);
		dE_factor_orig = (E_factor - old_E_factor);
		dE_factor = dE_factor_orig;

		if(sqrt(dE_factor_orig*dE_factor_orig) > (max_dE*TOTAL_ATOMS*k_b*T))
		dE_factor = (dE_factor_orig/sqrt(dE_factor_orig*dE_factor_orig))*max_dE*TOTAL_ATOMS*k_b*T;

		if((debug) && (mpi_id == 0))
		std::cout << " WARNING ! Probability factor " << P_factor << " greater than 1 ! " << " at step " << step
			  << " Actual pobability: " << E_Prob[E_id] << " Reference probability: " << E_Prob_ref[E_id] << std::endl;	
                }
		}
		else {
		if((debug) && (mpi_id == 0))
		std::cout << " WARNING ! Energy " << E 
		          << " is outside the range " 
		          << E_bins[0] << " - " << E_bins[bins-1] << " at step " << step << std::endl;
		}
		}

		old_E_factor = E_factor;

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
				double dx_mod = sqrt(dx*dx);

//				if((step > 0) && (step % tau == 0)){
				if((dx_mod > 0) && (dr_mod > 0))
				force[id][d] += (-dE_factor/(double(TOTAL_ATOMS)*dr_mod))*(dx/dr_mod);
				else{	
				if((debug) && (mpi_id == 0))
				std::cout << " WARNING ! particle did not move ! " 
					  << " dim = " << d << " dx_mod = " <<  dx_mod << " dr_mod " << dr_mod << std::endl;
				throw std::runtime_error(std::string("CFES simulation failed!"));
				}
//				}
				old_xyz(id, d) = xyz[id][d];
			}
		}

		//==================================================
		//Output
		//==================================================

		if(mpi_id == 0){
		if(step % output_steps == 0){
		e_factor_file << step << "  " << E << "  " << P_factor << "  " 
			      << E_factor << "  " << E_factor/(k_b*T) << "  " 
			      << dE_factor_orig << "  " << dE_factor << std::endl;

//		e_factor_bath_file << step << "  " << E_bath << "  " << P_factor_bath << "  " 
//				   << E_factor_bath << "  " << E_factor_bath/(k_b*T) << "  " 
				   
		e_factor_bath_file << step << "  " << C_Q_ref << "  " << C_Q << "  " << zeta_S << std::endl;

		for(int j = 0; j < bins; j++){
			probability_file << E_bins[j] << "  " << E_Prob[j] << std::endl;
//			probability_bath_file << E_bins[j] << "  " << E_Prob_bath[j] << std::endl;
		}	
		probability_file << "\n" << std::endl;
//		probability_bath_file << "\n" << std::endl;
		}}
		//==================================================
		}
	
	lammps_input_file_address = getFileAddress(directory, std::string("in.lammps_cfes.unfix"));
	lmp->input->file(lammps_input_file_address.c_str());

//==================================================

	delete lmp;
	MPI_Finalize();
}
