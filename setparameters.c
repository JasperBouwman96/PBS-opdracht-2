#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

void set_parameters(struct Parameters *p_parameters)
/* Set the parameters of this simulation */
{
  p_parameters->N        = 2000;          //number of particles
  p_parameters->numsteps = 20000;         //number of time steps
  p_parameters->kT       = 1.0;           //thermal energy
  p_parameters->mass     = 1.0;           //mass of a particle
  p_parameters->epsilon  = 1.0;           //LJ interaction strength
  p_parameters->sigma    = 1.0;           //LJ particle diameter
  p_parameters->rcut     = 2.5;           //cut-off distance for LJ interaction
  p_parameters->dt       = 0.001;         //integration time step
  p_parameters->L        = (struct Vec3D) {14.938, 14.938, 14.938};    //box size
  p_parameters->tau      = 0.5;           //Berendsen thermostat relaxation time
  p_parameters->rshell   = 0.4;           //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 500;
  strcpy(p_parameters->pdb_filename, "trajectories_short_1.0"); //filename (without extension) for pdb file
  p_parameters->num_dt_gr = 200;          //number of time steps between g_of_r updates
  p_parameters->num_dt_gr_start = 10000;  //number of time steps after which g_of_r updates start
  p_parameters->rmax_gr = 6.0;            //maximum r for g(r)
  p_parameters->num_bin_gr = 500;         //number of bins for g(r)
  strcpy(p_parameters->gr_filename, "gr_short_1.0.dat");  //filename for file containing the computed radial distribution function
  p_parameters->num_dt_msd = 1;           //number of time steps between msd updates
  p_parameters->num_dt_msd_start = 10000; //number of time steps after which msd updates start
  p_parameters->num_bin_msd = 500;        //number of bins for msd
  strcpy(p_parameters->msd_filename, "msd_short_1.0.dat"); //filename for file containing the computed msd

  if (p_parameters->rcut > p_parameters->L.x/2.0) fprintf(stderr,"Warning! rcut > Lx/2");
  if (p_parameters->rcut > p_parameters->L.y/2.0) fprintf(stderr,"Warning! rcut > Ly/2");
  if (p_parameters->rcut > p_parameters->L.z/2.0) fprintf(stderr,"Warning! rcut > Lz/2");
}
