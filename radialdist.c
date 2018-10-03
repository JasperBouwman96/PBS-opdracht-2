#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include <math.h>
#include "constants.h"

void radialdist(struct Nbrlist *p_nbrlist, struct Parameters *p_parameters)
{
    int k;
    int i;
    double rij;
    int bin;
    unsigned int N_bin = p_parameters->num_bin_gr;
    unsigned int gr[N_bin];
    double norm_gr[N_bin];
    double binwidth=0.5*p_parameters->L.x/N_bin;
    int count = 0;
    double surface;

for (i=0; i < N_bin; i++)
    {
       gr[i] = 0;
    }

    for(k = 0; k < p_parameters->N; ++k)
        {
          rij = sqrt(p_nbrlist->nbr[k].rij.sq);

         bin=rij/binwidth;
         gr[bin]++;
         count++;
    }
    int rho= (p_parameters->N/(p_parameters->L.x*p_parameters->L.y*p_parameters->L.z));

    FILE *F;
       F =  fopen("radial_dist.txt", "w");

   for (i=0;i<N_bin-1;i++)
    {
        surface=(4/3*PI*(((i+1.0)*binwidth)*((i+1.0)*binwidth))*((i+1.0)*binwidth)-(i*i*i))*binwidth*binwidth*binwidth;
        norm_gr[i]=gr[i]/(rho*surface*count*p_parameters->N);
        fprintf(F,"%f\n", norm_gr[i]);
    }



        fclose(F);
}






