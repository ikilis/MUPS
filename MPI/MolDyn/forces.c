#include "mpi.h"
#include <stdio.h>
#define MASTER 0
/*
 *  Compute forces and accumulate the virial and the potential
 */
extern double epot, vir;

void forces(int npart, double x[], double f[], double side, double rcoff)
{
  int i, j;
  double sideh, rcoffs;
  double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
  double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
  double forcex, forcey, forcez;

  vir = 0.0;
  epot = 0.0;
  sideh = 0.5 * side;
  rcoffs = rcoff * rcoff;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int start, end, chunk;
  chunk = (npart*3+size-1)/size;
  start = rank*chunk;
  end = (start+chunk)>(npart*3)?(npart*3):(start+chunk);

  MPI_Bcast(x, npart*3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(f, npart*3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);


  // printf("rank = %d, size = %d, start = %d, end = %d, chunk = %d \n", rank, size, start, end, chunk);

    // for (i = start; i < end; i += 3)
  for (i = 3*rank; i < npart*3; i += 3*size)
  {
    xi = x[i];
    yi = x[i + 1];
    zi = x[i + 2];
    fxi = 0.0;
    fyi = 0.0;
    fzi = 0.0;

    for (j = i + 3; j < npart * 3; j += 3)
    {
      xx = xi - x[j];
      yy = yi - x[j + 1];
      zz = zi - x[j + 2];
      if (xx < -sideh)
        xx += side;
      if (xx > sideh)
        xx -= side;
      if (yy < -sideh)
        yy += side;
      if (yy > sideh)
        yy -= side;
      if (zz < -sideh)
        zz += side;
      if (zz > sideh)
        zz -= side;
      rd = xx * xx + yy * yy + zz * zz;

      if (rd <= rcoffs)
      {
        rrd = 1.0 / rd;
        rrd2 = rrd * rrd;
        rrd3 = rrd2 * rrd;
        rrd4 = rrd2 * rrd2;
        rrd6 = rrd2 * rrd4;
        rrd7 = rrd6 * rrd;
        epot += (rrd6 - rrd3);
        r148 = rrd7 - 0.5 * rrd4;
        vir -= rd * r148;
        forcex = xx * r148;
        fxi += forcex;
        f[j] -= forcex;
        forcey = yy * r148;
        fyi += forcey;
        f[j + 1] -= forcey;
        forcez = zz * r148;
        fzi += forcez;
        f[j + 2] -= forcez;
      }
    }
    f[i] += fxi;
    f[i + 1] += fyi;
    f[i + 2] += fzi;
  }

  double reduced_f[npart*3], reduced_vir, reduced_epot;


  MPI_Reduce(f, reduced_f, npart*3, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&vir, &reduced_vir, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&epot, &reduced_epot, 1, MPI_DOUBLE, MPI_SUM,MASTER, MPI_COMM_WORLD);

  // MPI_Allreduce(f, reduced_f, npart*3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce(&vir, &reduced_vir, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce(&epot, &reduced_epot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  
  if(rank == MASTER)
  {
    epot = reduced_epot;
    vir = reduced_vir;

    for(int i = 0; i<npart*3; i++)
        f[i] = reduced_f[i];
  }

}
