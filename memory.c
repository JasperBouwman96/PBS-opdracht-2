#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "nbrlist.h"

void alloc_vectors(struct Vectors *p_vectors, size_t N)
/* Allocate the arrays in 'vectors' needed to store information of all particles */
{
    p_vectors->ID  = (size_t *) malloc(N*sizeof(size_t));
    p_vectors->r  = (struct Vec3D *) malloc(N*sizeof(struct Vec3D));
    p_vectors->dr = (struct Vec3D *) malloc(N*sizeof(struct Vec3D));
    p_vectors->v  = (struct Vec3D *) malloc(N*sizeof(struct Vec3D));
    p_vectors->f  = (struct Vec3D *) malloc(N*sizeof(struct Vec3D));
}

void free_vectors(struct Vectors *p_vectors)
/* Free the arrays in 'vectors' */
{
    free(p_vectors->ID);
    p_vectors->ID = NULL;
    free(p_vectors->r);
    p_vectors->r = NULL;
    free(p_vectors->dr);
    p_vectors->dr = NULL;
    free(p_vectors->v);
    p_vectors->v = NULL;
    free(p_vectors->f);
    p_vectors->f = NULL;
}

void alloc_memory(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Allocate all variables needed in the MD simulation */
{
    alloc_vectors(p_vectors, p_parameters->N);
    alloc_nbrlist(p_parameters, p_nbrlist);
}

void free_memory(struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Free the memory allocated by alloc_memory */
{
    free_vectors(p_vectors);
    free_nbrlist(p_nbrlist);
}
