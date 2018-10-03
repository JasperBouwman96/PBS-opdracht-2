#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist)
/* Allocate arrays needed to store the cell-linked-list data*/
{
    size_t Mx, My, Mz, Mtot;
    const double rlist = p_parameters->rcut+p_parameters->rshell;
    const size_t N = p_parameters->N;
    Mx = floor(p_parameters->L.x/rlist);
    My = floor(p_parameters->L.y/rlist);
    Mz = floor(p_parameters->L.z/rlist);
    p_celllist->Mx = Mx;
    p_celllist->My = My;
    p_celllist->Mz = Mz;
    Mtot = Mx*My*Mz;
    p_celllist->head = (size_t *) malloc(Mtot*sizeof(size_t));
    p_celllist->num_cells_max = Mtot;
    p_celllist->num_cells = Mtot;
    p_celllist->particle2cell = (size_t *) malloc(N*sizeof(size_t));
    p_celllist->list = (size_t *) malloc(N*sizeof(size_t));
    p_celllist->num_cells_max = N;
}

void free_celllist(struct Celllist *p_celllist)
/* Free arrays used for the cell-linked-list */
{
    free(p_celllist->head);
    p_celllist->head = NULL;
    free(p_celllist->particle2cell);
    p_celllist->particle2cell = NULL;
    free(p_celllist->list);
    p_celllist->list = NULL;
}

void build_celllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Celllist *p_celllist)
/* Build the cell-linked-list */
{
    size_t Mx, My, Mz, Mtot;
    size_t i, ix, iy, iz, icell;
    const double rlist = p_parameters->rcut+p_parameters->rshell;
    struct Vec3D mL;
    const size_t N = p_parameters->N;
    size_t *particle2cell, *head, *celllist;
    struct Vec3D *r;

    Mx = floor(p_parameters->L.x/rlist);
    My = floor(p_parameters->L.y/rlist);
    Mz = floor(p_parameters->L.z/rlist);
    p_celllist->Mx = Mx;
    p_celllist->My = My;
    p_celllist->Mz = Mz;
    Mtot = Mx*My*Mz;
    if (Mtot > p_celllist->num_cells_max)
    {
        p_celllist->head = (size_t *) realloc(p_celllist->head, Mtot*sizeof(size_t));
        p_celllist->num_cells_max = Mtot;
    }
    if (N > p_celllist->num_part_max)
    {
        p_celllist->particle2cell = (size_t *) realloc(p_celllist->particle2cell, N*sizeof(size_t));
        p_celllist->list = (size_t *) realloc(p_celllist->list, N*sizeof(size_t));
        p_celllist->num_part_max = N;
    }
    mL.x = ((double) Mx)/p_parameters->L.x;
    mL.y = ((double) My)/p_parameters->L.y;
    mL.z = ((double) Mz)/p_parameters->L.z;
    head = p_celllist->head;
    for(icell=0; icell < Mtot; ++icell)
        head[icell] = SIZE_MAX;
    particle2cell = p_celllist->particle2cell;
    celllist = p_celllist->list;
    r = p_vectors->r;
    for(i=0; i < N; ++i)
    {
        ix = floor(r[i].x*mL.x);
        iy = floor(r[i].y*mL.y);
        iz = floor(r[i].z*mL.z);
        icell = ix + Mx * (iy + iz*My);
        particle2cell[i] = icell;
        celllist[i] = head[icell];
        head[icell] = i;
    }
}

void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist)
/* Allocate arrays needed to store the neighbor list */
{
    double Ndouble = (double) p_parameters->N;
    double rlist = p_parameters->rcut + p_parameters->rshell;
    // Estimate the number of neighbors using the average density
    double Nnbr = 0.5*Ndouble*(4.0/3.0*PI*rlist*rlist*rlist) / (p_parameters->L.x*p_parameters->L.y*p_parameters->L.z);
    size_t num_nbrs_max = ceil(0.6*Ndouble*(Nnbr + 2.0*sqrt(Nnbr)));
    p_nbrlist->p_celllist = (struct Celllist *) malloc(sizeof(struct Celllist));
    alloc_celllist(p_parameters, p_nbrlist->p_celllist);
    p_nbrlist->num_nbrs_max = num_nbrs_max;
    p_nbrlist->nbr = (struct Pair *) malloc(num_nbrs_max*sizeof(struct Pair));
    p_nbrlist->dr = (struct DeltaR *) malloc(p_parameters->N*sizeof(struct DeltaR));
}

void free_nbrlist(struct Nbrlist *p_nbrlist)
/* Free arrays use to store the neighbor list */
{
    free_celllist(p_nbrlist->p_celllist);
    free(p_nbrlist->p_celllist);
    p_nbrlist->p_celllist = NULL;
    free(p_nbrlist->nbr);
    p_nbrlist->nbr = NULL;
    free(p_nbrlist->dr);
    p_nbrlist->dr = NULL;
}

void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Build the neighbor list */
{
    size_t Mx, My, Mz;
    size_t ix, iy, iz, icell;
    size_t ix_nbr, iy_nbr, iz_nbr, inbr;
    size_t num_nbrs, num_nbrs_max;
    size_t i, j, k;
    const double rlist = p_parameters->rcut+p_parameters->rshell; /* the radius for inclusion in the list is rcut + rshell */
    const double rlist_sq = rlist*rlist;
    struct Vec3D ri;
    struct Vec3D *r = p_vectors->r;
    struct DeltaR rij;
    const struct DeltaR dr = {0.0,0.0,0.0,0.0};
    size_t *head, *particle2cell, *celllist;
    const int nbr_indcs[13][3] = {{0,0,1},{0,1,-1},{0,1,0},{0,1,1},{1,-1,-1},{1,-1,0},{1,-1,1},{1,0,-1},
        {1,0,0},{1,0,1},{1,1,-1},{1,1,0},{1,1,1}
    };
    size_t N = p_parameters->N;
    struct Pair *nbr;

    // First build a cell-linked-list
    build_celllist(p_parameters, p_vectors, p_nbrlist->p_celllist);

    /*  Use the cell-linked-list to build a neighbor list.
        Pairs are included to the neighbor list of their distance is less than rcut+rshell. */
    Mx = p_nbrlist->p_celllist->Mx;
    My = p_nbrlist->p_celllist->My;
    Mz = p_nbrlist->p_celllist->Mz;
    num_nbrs_max = p_nbrlist->num_nbrs_max;
    num_nbrs = 0;
    nbr = p_nbrlist->nbr;
    particle2cell = p_nbrlist->p_celllist->particle2cell;
    head = p_nbrlist->p_celllist->head;
    celllist = p_nbrlist->p_celllist->list;
    for(i = 0 ; i < N ; ++i)
    {
        // find neigbors of particle i in its own cell
        ri = r[i];
        for( j = celllist[i]; j != SIZE_MAX; j = celllist[j])
        {
            rij.x = ri.x - r[j].x;
            rij.y = ri.y - r[j].y;
            rij.z = ri.z - r[j].z;
            rij.sq = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
            if(rij.sq < rlist_sq)
            {
                if (num_nbrs >= num_nbrs_max)
                {
                    num_nbrs_max += 5*p_parameters->N;
                    nbr = (struct Pair *) realloc(nbr, num_nbrs_max*sizeof(struct Pair));
                    p_nbrlist->nbr = nbr;
                    p_nbrlist->num_nbrs_max = num_nbrs_max;
                }
                nbr[num_nbrs].i = i;
                nbr[num_nbrs].j = j;
                nbr[num_nbrs].rij = rij;
                num_nbrs++;
            }
        }

        // next find neighbors of particle i in 1 of its 13 neighboring cells
        icell = particle2cell[i];
        ix = icell % Mx;
        icell = icell / Mx;
        iy = icell % My;
        iz = icell / My;
        for( k=0; k< 13; ++k)
        {
            ri = p_vectors->r[i];
            ix_nbr = (ix + nbr_indcs[k][0]);
            iy_nbr = (iy + nbr_indcs[k][1]);
            iz_nbr = (iz + nbr_indcs[k][2]);
            // The if-statements below implement periodic boundary conditions
            if (ix_nbr == SIZE_MAX)
            {
                ri.x += p_parameters->L.x;
                ix_nbr = Mx-1;
            }
            else if (ix_nbr>=Mx)
            {
                ri.x -= p_parameters->L.x;
                ix_nbr = 0;
            }
            if (iy_nbr == SIZE_MAX)
            {
                ri.y += p_parameters->L.y;
                iy_nbr = My-1;
            }
            else if (iy_nbr>=My)
            {
                ri.y -= p_parameters->L.y;
                iy_nbr = 0;
            }
            if (iz_nbr == SIZE_MAX)
            {
                ri.z += p_parameters->L.z;
                iz_nbr = Mz-1;
            }
            else if (iz_nbr>=Mz)
            {
                ri.z -= p_parameters->L.z;
                iz_nbr = 0;
            }
            inbr = ix_nbr + Mx * (iy_nbr + iz_nbr*My);
            for(j = head[inbr]; j != SIZE_MAX; j = celllist[j])
            {
                rij.x = ri.x - p_vectors->r[j].x;
                rij.y = ri.y - p_vectors->r[j].y;
                rij.z = ri.z - p_vectors->r[j].z;
                rij.sq = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
                if(rij.sq < rlist_sq)
                {
                    if (num_nbrs >= num_nbrs_max)
                    {
                        num_nbrs_max += 5*p_parameters->N;
                        nbr = (struct Pair *) realloc(nbr, num_nbrs_max*sizeof(struct Pair));
                        p_nbrlist->nbr = nbr;
                        p_nbrlist->num_nbrs_max = num_nbrs_max;
                    }
                    nbr[num_nbrs].i = i;
                    nbr[num_nbrs].j = j;
                    nbr[num_nbrs].rij = rij;
                    num_nbrs++;
                }
            }
        }
    }
    p_nbrlist->num_nbrs = num_nbrs;
    p_nbrlist->dr = (struct DeltaR *) realloc(p_nbrlist->dr, N * sizeof(struct DeltaR));
    for(i=0; i < N; ++i) /*initialize particle displacements (with respect to creation time) to zero */
        p_nbrlist->dr[i] = dr;
}

int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Update the connecting vectors of the neighbor list and if needed rebuild it.*/
{
    size_t i,j,k;
    const double dr_sq_max = 0.25*(p_parameters->rshell*p_parameters->rshell);
    int isRebuild = 0;
    struct DeltaR rij;
    struct Pair * nbr = p_nbrlist->nbr;
    struct Vec3D * dr = p_vectors->dr;
    // The neighbor list needs to be rebuild if one of the particles has displaced more then 0.5*rshell
    for(i=0; i < p_parameters->N; ++i)
        if( (p_nbrlist->dr[i].sq) > dr_sq_max)
        {
            isRebuild = 1;
        }
    if (isRebuild) // rebuild neighbor list
        build_nbrlist(p_parameters, p_vectors, p_nbrlist);
    else // If no rebuild is needed, update the values of the connecting vectors
    {
        for(k = 0; k < p_nbrlist->num_nbrs; ++k)
        {
            i = nbr[k].i;
            j = nbr[k].j;
            rij = nbr[k].rij;
            rij.x += (dr[i].x-dr[j].x);
            rij.y += (dr[i].y-dr[j].y);
            rij.z += (dr[i].z-dr[j].z);
            rij.sq = rij.x*rij.x+rij.y*rij.y+rij.z*rij.z;
            p_nbrlist->nbr[k].rij = rij;
        }
    }
    return isRebuild;
}
