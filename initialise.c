#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

void initialise_variables(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t *p_step, double *p_time)
/* Initialize all the variables of an MD simulation. */
{
    size_t i;
    srand(13); // Positive integer as seed for random number generator
    for (i=0; i<p_parameters->N; i++)
        p_vectors->ID[i] = i; // Give each particle a unique ID. Note this ID is actually not used in this code.
    initialise_positions(p_parameters, p_vectors);
    initialise_velocities(p_parameters, p_vectors);
    *p_step = 0; // Initialize the step as zero
    *p_time = 0.0; // Initialize the time to zero
    return;
}

void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors)
/*  Initialize positions of all particles.
    Particles are initialized on a square lattice.
*/
{
    struct Vec3D dr;
    double dl;
    int nx,ny,nz,i,j,k,ipart;

    dl = pow(p_parameters->L.x*p_parameters->L.y*p_parameters->L.z/((double) p_parameters->N),1.0/3.0);
    nx = (int) ceil(p_parameters->L.x/dl);
    ny = (int) ceil(p_parameters->L.y/dl);
    nz = (int) ceil(p_parameters->L.z/dl);
    dr.x = p_parameters->L.x/(double)nx;
    dr.y = p_parameters->L.y/(double)ny;
    dr.z = p_parameters->L.z/(double)nz;
    ipart = 0;
    for(i=0; i<nx; ++i)
        for(j=0; j<ny; ++j)
            for(k=0; k<nz; ++k, ++ipart)
            {
                if (ipart >= p_parameters->N) break;
                p_vectors->r[ipart].x = (i+0.5)*dr.x;
                p_vectors->r[ipart].y = (j+0.5)*dr.y;
                p_vectors->r[ipart].z = (k+0.5)*dr.z;
//      p_vectors->r[ipart].x = p_parameters->L.x*generate_uniform_random();
//      p_vectors->r[ipart].y = p_parameters->L.y*generate_uniform_random();
//      p_vectors->r[ipart].z = p_parameters->L.z*generate_uniform_random();

            }
}

void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
/*  Initialize the velocities of all particles.
    Initial velocities are sampled according to the Maxwell-Boltzmann distribution
*/
{
    double sqrtktm;
    struct Vec3D sumv;
    int i;

    sqrtktm = sqrt(p_parameters->kT/p_parameters->mass);
    sumv   = (struct Vec3D)
    {
        0.0
    };

    for (i=0; i<p_parameters->N; i++)
    {
        p_vectors->v[i].x = sqrtktm*gauss();
        p_vectors->v[i].y = sqrtktm*gauss();
        p_vectors->v[i].z = sqrtktm*gauss();
        sumv.x += p_vectors->v[i].x;
        sumv.y += p_vectors->v[i].y;
        sumv.z += p_vectors->v[i].z;
    }

    sumv.x /= ((double)(p_parameters->N)); /* remove average velocity */
    sumv.y /= ((double)(p_parameters->N)); /* so total momentum is zero */
    sumv.z /= ((double)(p_parameters->N));
    for (i=0; i<p_parameters->N; i++)
    {
        p_vectors->v[i].x -= sumv.x;
        p_vectors->v[i].y -= sumv.y;
        p_vectors->v[i].z -= sumv.z;
    }
}
