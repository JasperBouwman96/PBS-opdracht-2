#ifndef TYPES_MD_H_
#define TYPES_MD_H_

/* This header file contains definitions of struct types used in the molecular dynamics code */

struct Vec3D
/* Structure to store x, y, and z component of a 3D vector */
{
    double x, y, z;
};

struct Parameters
/* Structure to store all parameters. These parameters are set by the function set_parameters. */
{
    size_t N;           /* Number of particles */
    size_t numsteps;    /* Number of time steps */
    double  kT;         /* Thermal energy */
    double  mass;       /* Mass of a particle */
    double  epsilon;    /* LJ interaction strength */
    double  sigma;      /* LJ particle diameter */
    double  rcut;       /* Cut-off distance for LJ interaction */
    double  dt;         /* integration time step */
    struct  Vec3D L;    /* Box sizes in 3 direction*/
    double  tau;        /* Berendsen thermostat relaxation time */
    double  rshell;     /* Shell thickness for neighbor list */
    size_t  num_dt_pdb; /* Number of time steps between pdb saves */
    char    pdb_filename[1024]; /* filename (without extension) for pdb file */
    size_t  num_dt_gr;  /* Number of time steps between g_of_r updates */
    size_t  num_dt_gr_start;    /* Number of time steps after which g_of_r updates start */
    double  rmax_gr;    /* maximum r for g(r) */
    size_t  num_bin_gr; /* number of bins for g(r) */
    char    gr_filename[1024];  /* filename to save g_of_r data */
    size_t  num_dt_msd;  /* Number of time steps between msd updates */
    size_t  num_dt_msd_start;   /* Number of time steps after which msd updates start */
    size_t  num_bin_msd;         /* number of bins for msd */
    char    msd_filename[1024]; /* filename to save msd data */
};

struct DeltaR
/* Structure to store a 3D vector and its square length. */
{
    double x, y, z;     /* x, y and z coordinates */
    double sq;          /* square length */
};

struct Vectors
/* Structure with pointers to all arrays relevant for a MD simulation */
{
    size_t *ID;              /* unique ID */
    struct Vec3D *r;         /* positions */
    struct Vec3D *dr;        /* displacements */
    struct Vec3D *v;         /* velocities */
    struct Vec3D *f;         /* forces */
};

struct Pair
/* Structure to store a pair of particles: its indices and connecting vector*/
{
    size_t i,j;
    struct DeltaR rij;
};

struct Celllist
/* Structure to store a cell-linked-list */
{
    size_t *head;
    size_t *list;
    size_t *particle2cell;
    size_t num_cells, num_cells_max, num_part_max;
    size_t Mx, My, Mz;
};

struct Nbrlist
/* Structure to store a nbr-list */
{
    struct Celllist * p_celllist;
    size_t num_nbrs, num_nbrs_max;
    struct Pair *nbr;
    struct DeltaR *dr;   /* displacements particles with respect to nbrlist creation time */
};

#endif /* TYPES_MD_H_ */
