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

  MPI_Bcast(x, npart*3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(f, npart*3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

  // printf("rank = %d, size = %d \n", rank, size);

  MPI_Status status;

  if(rank == MASTER)
  {

    // posalji im instrukcije
    int pomeraj = 3*(size-1);
    for(int i =1; i < size; i++){
      // printf("MASTER salje pomeraj\n");
      MPI_Send(&pomeraj, 1, MPI_INT, i, 103, MPI_COMM_WORLD);
    }
    double reduced_f[npart*3], reduced_epot, reduced_vir;

    // primi podatke nazad
    for(int i=1; i<size;i++){
      // printf("MASTER CEKA ODGOVOR f \n");
      MPI_Recv(reduced_f, npart*3, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
      // printf("MASTER CEKA ODGOVOR epot \n");
      MPI_Recv(&reduced_epot, 1, MPI_DOUBLE, i, 101, MPI_COMM_WORLD, &status);
      // printf("MASTER CEKA ODGOVOR vir \n");
      MPI_Recv(&reduced_vir, 1, MPI_DOUBLE, i, 102, MPI_COMM_WORLD, &status);

      // printf("MASTER primio sve odgovore od %d \n", i);
      epot+=reduced_epot;
      vir+=reduced_vir;
      for(int j=0;j<npart*3;j++)f[j]+=reduced_f[j];
      // printf("MASTER azurirao matricu od %d rob \n", i);
    }

  }else
  {
    // primi pomeraj
    int pomeraj;
    // printf("%d CEKA pomeraj \n", rank);
    MPI_Recv(&pomeraj, 1, MPI_INT, MASTER, 103, MPI_COMM_WORLD, &status);

    int tmpRank = rank-1;

    for (i = 3*tmpRank; i < npart * 3; i += pomeraj)
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
    
    // posalji nazad
    // printf("%d salje sve masteru \n", rank);
    MPI_Send(f, npart*3, MPI_DOUBLE, MASTER, 100, MPI_COMM_WORLD);
    MPI_Send(&epot, 1, MPI_DOUBLE, MASTER, 101, MPI_COMM_WORLD);
    MPI_Send(&vir, 1, MPI_DOUBLE, MASTER, 102, MPI_COMM_WORLD);
  }

}
