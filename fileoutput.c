#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "structs.h"

void record_trajectories(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
/*  Write the particle positions to a pdf file
    The filename (without extension) is given by p_parameters->pdb_filename.
    If reset = 1 the data is written to the file deleting data it possibly contained.
    If reset = 0 the data is appended. */
{
  size_t i;
  FILE* fp_traj;
  char filename[1024];

  snprintf(filename, 1024, "%s%s", p_parameters->pdb_filename, ".pdb");
  if (reset == 1)
    {
      fp_traj = fopen( filename, "w" );
    }
  else
    {
      fp_traj = fopen( filename, "a" );
    }

  fprintf( fp_traj, "REMARK TIME = %f\n",time);
  fprintf( fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n",p_parameters->L.x,p_parameters->L.y,p_parameters->L.z,90.0,90.0,90.0,"P 1","1");
  for ( i=0 ; i<p_parameters->N ; i++ )
    {
        fprintf( fp_traj, "HETATM %5u  C 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", (unsigned int) i%100000, p_vectors->r[i].x, p_vectors->r[i].y, p_vectors->r[i].z);
    }
  fprintf( fp_traj, "END\n");

  fclose(fp_traj);
}
