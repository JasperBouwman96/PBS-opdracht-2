#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
/* Compute forces on particles using the pairs in a neighbor list.
This function returns the total potential energy of the system. */
{
    size_t i,j,k;
    struct Vec3D df;
    double rcutsq,sigmasq,sr2,sr6,sr12,fr,prefctr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t N = p_parameters->N;

    rcutsq = p_parameters->rcut*p_parameters->rcut;
    sigmasq = p_parameters->sigma*p_parameters->sigma;

    double Epot = 0.0, Epot_cutoff;
    sr2 = sigmasq/rcutsq;
    sr6 = sr2*sr2*sr2;
    sr12 = sr6*sr6;
    Epot_cutoff = sr12 - sr6;

    for (i=0; i < N; i++)
    // initialize the forces to zero
        f[i] = (struct Vec3D) {0.0, 0.0, 0.0}; /*initialize forces to zero*/

    for (k=0; k < num_nbrs; k++)
    {
        // for each pair in the neighbor list compute the pair forces
        rij = nbr[k].rij;
        i = nbr[k].i;
        j = nbr[k].j;
        if (rij.sq < rcutsq)
        // Compute forces if the distance is smaller than the cutoff distance
        {
            // pair forces are given by the LJ interaction
            sr2 = sigmasq/rij.sq;
            sr6   = sr2*sr2*sr2;
            sr12  = sr6*sr6;
            Epot += (sr12 - sr6 - Epot_cutoff);
            fr    = (2.0*sr12-sr6)/rij.sq; //force divided by distance
            df.x = fr*rij.x;
            df.y = fr*rij.y;
            df.z = fr*rij.z;
            f[i].x += df.x;
            f[i].y += df.y;
            f[i].z += df.z;
            f[j].x -= df.x;
            f[j].y -= df.y;
            f[j].z -= df.z;
        }
    }

    prefctr = 24.0*p_parameters->epsilon;
    for (i=0; i<p_parameters->N; i++)
    //multiply with prefactors
    {
        f[i].x *= prefctr;
        f[i].y *= prefctr;
        f[i].z *= prefctr;
    }
    Epot *= 4.0*p_parameters->epsilon;

    return Epot;
}
