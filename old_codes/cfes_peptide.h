#include "includeallglobals.h"

void lammps_equil(){
	lammps_commands_string(LAMMPS_HANDLE,
	"compute	E peptide pe/atom pair bond angle dihedral improper kspace \n"
	"compute	Etot peptide reduce sum c_E \n"
	"compute 	coord peptide property/atom x y z \n"
	"compute        force peptide property/atom fx fy fz \n"
	"variable	pepE equal c_Etot \n"
	"dump   	traj peptide custom ${writestep} ${sname}_cfes_md_equil.lammpstrj id type xu yu zu \n" 
	"dump_modify	traj append yes \n"
        "fix    	file all print ${writestep} '$s $T $P ${kin} ${pot} ${pepE}' append 'thermo_cfes_md_equil.dat' screen 'no' title '#' \n"
	"fix 		hist all print 1 '${pepE}' 'hist.dat' screen 'no' title '#' \n"
	"fix		constraint peptide shake 0.0001 20 0 m 1.0080 \n"
	"fix		pepmd peptide npt temp ${TrgtTemp} ${TrgtTemp} 100.0 iso ${TrgtPress} ${TrgtPress} 100.0 \n"
	"fix    	mom all momentum 100 linear 1 1 1 angular rescale \n"
	"run		1 \n");
}

void lammps_cfes(){
	lammps_commands_string(LAMMPS_HANDLE,
	"compute	E peptide pe/atom pair bond angle dihedral improper kspace \n"
	"compute	Etot peptide reduce sum c_E \n"
	"compute 	coord peptide property/atom x y z \n"
	"compute        force peptide property/atom fx fy fz \n"
	"variable	pepE equal c_Etot \n"
	"dump   	traj peptide custom ${writestep} ${sname}_cfes_md.lammpstrj id type xu yu zu \n" 
	"dump_modify	traj append yes \n"
	"dump           trajF peptide custom ${writestep} ${sname}_cfes_md_forces.lammpstrj id type fx fy fz \n"
	"dump_modify	trajF append yes \n"
	 "fix    	file all print ${writestep} '$s $T $P ${kin} ${pot} ${pepE}' append 'thermo_cfes_md_equil.dat' screen 'no' title '#' \n"
	"fix 		hist all print 1 '${pepE}' 'hist.dat' screen 'no' title '#' \n"
	"fix 		xcoord all print 1 'c_coord[1]' 'xcoord.dat' screen 'no' title '#' \n"
	"fix 		ycoord all print 1 'c_coord[2]' 'ycoord.dat' screen 'no' title '#' \n"
	"fix 		zcoord all print 1 'c_coord[3]' 'zcoord.dat' screen 'no' title '#' \n"
	"fix		constraint peptide shake 0.0001 20 0 m 1.0080 \n"
	"fix		pepmd peptide npt temp ${TrgtTemp} ${TrgtTemp} 100.0 iso ${TrgtPress} ${TrgtPress} 100.0 \n"
	"fix    	mom all momentum 100 linear 1 1 1 angular rescale \n"
	"run		1 \n");
}

void lammps_cfes_addforce(){
	lammps_commands_string(LAMMPS_HANDLE,
	"variable 	newfx atomfile newfx.dat \n"
	"variable 	newfy atomfile newfy.dat \n"
	"variable 	newfz atomfile newfz.dat \n"
	"fix		qstat peptide addforce v_newfx v_newfy v_newfz \n"
	"run		1 \n"
	"unfix  	qstat \n");
}

void lammps_unfix(){
	lammps_commands_string(LAMMPS_HANDLE,
	"unfix	pepmd \n"
	"unfix	mom \n"
	"unfix	file \n"
	"unfix	hist \n"
	"unfix  xcoord \n"
	"unfix  ycoord \n"
	"unfix  zcoord \n"
	"unfix	constraint \n"
	"undump	traj \n"
	"undump	trajF \n");
}
