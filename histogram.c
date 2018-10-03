#include <stdio.h>
#include <stdlib.h>
#include "structs.h"


void histogram(struct Vectors *p_vectors, struct Parameters *p_parameters)
{
    int i;
    double v_max = 4 , v_min = -4;
    int bin;
    unsigned int N_bin = 100;
    unsigned int Histo[N_bin];
    double binwidth=(v_max-v_min)/(double)N_bin;
    int count = 0;
    for (i=0; i < N_bin; i++)
    {
       Histo[i] = 0;
    }
    for (i=0; i < p_parameters->N; i++)
    {
        //number of bins =100 with a width of 8/100 = 0.08
        bin=(p_vectors->v[i].x - v_min)/binwidth;
        if (0<bin && bin<N_bin)
        {
           Histo[bin]++;
            count++;
        }

    }
    double norm_histo[N_bin];
    for (i=0; i<N_bin; i++)
    {
            norm_histo[i] = (double)Histo[i]/(binwidth*(double)count);
    }
   FILE *F;
   F=fopen("histogram.txt","w");
   for (i=0; i<N_bin; i++)
   {
       fprintf(F,"%12.8f\n", norm_histo[i]);
   }
   fclose(F);
}
//histogram(p_vectors v.x)


