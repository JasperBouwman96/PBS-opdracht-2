#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include <math.h>
#include "setparameters.h"
#include "constants.h"


int i,j;
int xmax = 2000;
int ymax = 3;
int count[3000] = {0};     //
double cor[3000] = {0.0};  //lenght ringbuffer
int ncor = 3000;

void initial(struct Vectors *p_vectors, struct Parameters *p_parameters,double *store, int *frame){
for (i=0; i<p_parameters->N; i++){
    store[(int)frame*xmax*ymax + 0*xmax + i] = p_vectors->r[i].x;
    store[(int)frame*xmax*ymax + 1*xmax + i] = p_vectors->r[i].y;
    store[(int)frame*xmax*ymax + 2*xmax + i] = p_vectors->r[i].z;
    frame++;
    printf("hoi");
}
}
void calc_msd(struct Vectors *p_vectors, struct Parameters *p_parameters,double *store, int *frame){
    int maxcor = ncor-1;
    if (frame < ncor){
        maxcor = frame;
    }
    int curframe = ((int)frame%ncor);
    int prevframe = ((curframe-1+ncor)%ncor);
    for (i=0; i<p_parameters->N; i++)
    {
    double deltax=p_vectors->r[i].x - store[prevframe*xmax*ymax + 0*xmax + i];
    double deltay=p_vectors->r[i].y - store[prevframe*xmax*ymax + 1*xmax + i];
    double deltaz=p_vectors->r[i].z - store[prevframe*xmax*ymax + 2*xmax + i];
    double xunfold = store[prevframe*xmax*ymax + 0*xmax + i]+deltax-(p_parameters->L.x*floor(deltax/p_parameters->L.x));
    double yunfold = store[prevframe*xmax*ymax + 1*xmax + i]+deltay-(p_parameters->L.y*floor(deltay/p_parameters->L.y));
    double zunfold = store[prevframe*xmax*ymax + 2*xmax + i]+deltaz-(p_parameters->L.z*floor(deltaz/p_parameters->L.z));
    store[curframe*xmax*ymax + 0*xmax + i] = xunfold;
    store[curframe*xmax*ymax + 0*xmax + i] = yunfold;
    store[curframe*xmax*ymax + 0*xmax + i] = zunfold;
    for (int icor = 0; icor<maxcor; icor++){
        int icorframe = ((curframe-icor+ncor)%ncor);
        cor[icor] = cor[icor] + (xunfold-store[icorframe*xmax*ymax + 0*xmax + i])*(xunfold-store[icorframe*xmax*ymax + 0*xmax + i])
        + (yunfold-store[icorframe*xmax*ymax + 1*xmax + i])*(yunfold-store[icorframe*xmax*ymax + 1*xmax + i])
        + (zunfold-store[icorframe*xmax*ymax + 2*xmax + i])*(zunfold-store[icorframe*xmax*ymax + 2*xmax + i]);
        count[icor]++;
    }
    }
}






