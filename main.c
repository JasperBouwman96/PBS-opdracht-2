/******************************************************************************/
/*                                                    			              */
/*  A Molecular Dynamics simulation of Lennard-Jones disks (2d)               */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*  Version 4.0, 18/9/2018                                                    */
/*                                                                            */
/*  Dr. Ir. J.T. Padding: version 1.1, 30/1/2013                              */
/*  Jeroen Hofman:        version 1.2, 28/7/2015                              */
/*                                                    			              */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"
#include "histogram.h"
#include "radialdist.h"
#include "msd.h"

int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    size_t step;
    double Ekin, Epot, time;

    set_parameters(&parameters);
    alloc_memory(&parameters,&vectors,&nbrlist);
    initialise_variables(&parameters,&vectors,&step,&time);

    build_nbrlist(&parameters, &vectors, &nbrlist);
    Epot = calculate_forces(&parameters,&nbrlist,&vectors);
    record_trajectories(1,&parameters,&vectors,time);

    //intitialise the matrix for msd
    double store[2000*3*3000] = {0};
    int frame = 0;

    while (step < parameters.numsteps)
    {
        step++;
        time += parameters.dt;

        Ekin = update_velocities_half_dt(&parameters,&nbrlist,&vectors);
        berendsen_thermostat(&parameters,&vectors,Ekin);
        update_positions(&parameters,&nbrlist,&vectors);
        Epot = calculate_forces(&parameters,&nbrlist,&vectors);
        boundary_conditions(&parameters,&vectors);
        update_nbrlist(&parameters, &vectors, &nbrlist);
        Ekin = update_velocities_half_dt(&parameters,&nbrlist,&vectors);

        //printf("Step %lu, Time %f, Epot %f, Ekin %f, Etot %f\n", (long unsigned) step, time, Epot, Ekin, Epot+Ekin);
        //if (step%parameters.num_dt_pdb ==0) record_trajectories(0,&parameters,&vectors,time);

        //make velocity histogram of v.x component
        /*if (step == 10000)
        {
            histogram(&vectors, &parameters);
        }*/

        // initialise first frame for msd
        if(step == 1599){
        initial(&vectors, &parameters, &store, &frame);
        }
    }
    free_memory(&vectors,&nbrlist);

    return 0;
}
