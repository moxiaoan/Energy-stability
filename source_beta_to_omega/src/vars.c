#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"
#include "ScLETD.h"
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif

void global_index(void)
{
  int i, j, k, gx, gy, gz;
  double x, y, z, cx, cy, cz;

  if (cart_id[0] == 0)
  {
    gx = 0;
  }
  else
  {
    gx = cart_id[0] * nx - 2 * nghost * cart_id[0];
  }
  if (cart_id[1] == 0)
  {
    gy = 0;
  }
  else
  {
    gy = cart_id[1] * ny - 2 * nghost * cart_id[1];
  }
  if (cart_id[2] == 0)
  {
    gz = 0;
  }
  else
  {
    gz = cart_id[2] * nz - 2 * nghost * cart_id[2];
  }
  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        fieldgx[k * nx * ny + j * nx + i] = gx + i;
        fieldgy[k * nx * ny + j * nx + i] = gy + j;
        fieldgz[k * nx * ny + j * nx + i] = gz + k;
      }
    }
  }
}

void init_para(void)
{
  if (left < 0)
  {
    ix1 = 0;
    ix2 = ix1;
    ix3 = ix1 + nx - nghost;
    ix4 = nx;
  }
  else if (right < 0)
  {
    ix1 = 0;
    ix2 = ix1 + nghost;
    ix3 = ix1 + nx;
    ix4 = nx;
  }
  else
  {
    ix1 = 0;
    ix2 = ix1 + nghost;
    ix3 = ix1 + nx - nghost;
    ix4 = nx;
  }
  lnx = ix4 - ix1;
  if (periodic)
  {
    gnx = lnx * procs[0] - 2 * nghost * procs[0];
  }
  else
  {
    gnx = lnx * procs[0] - 2 * nghost * (procs[0] - 1);
  }

  if (top < 0)
  {
    iy1 = 0;
    iy2 = iy1;
    iy3 = iy1 + ny - nghost;
    iy4 = ny;
  }
  else if (bottom < 0)
  {
    iy1 = 0;
    iy2 = iy1 + nghost;
    iy3 = iy1 + ny;
    iy4 = ny;
  }
  else
  {
    iy1 = 0;
    iy2 = iy1 + nghost;
    iy3 = iy1 + ny - nghost;
    iy4 = ny;
  }
  lny = iy4 - iy1;
  if (periodic)
  {
    gny = lny * procs[1] - 2 * nghost * procs[1];
  }
  else
  {
    gny = lny * procs[1] - 2 * nghost * (procs[1] - 1);
  }

  if (front < 0)
  {
    iz1 = 0;
    iz2 = iz1;
    iz3 = iz1 + nz - nghost;
    iz4 = nz;
  }
  else if (back < 0)
  {
    iz1 = 0;
    iz2 = iz1 + nghost;
    iz3 = iz1 + nz;
    iz4 = nz;
  }
  else
  {
    iz1 = 0;
    iz2 = iz1 + nghost;
    iz3 = iz1 + nz - nghost;
    iz4 = nz;
  }
  lnz = iz4 - iz1;
  if (periodic)
  {
    gnz = lnz * procs[2] - 2 * nghost * procs[2];
  }
  else
  {
    gnz = lnz * procs[2] - 2 * nghost * (procs[2] - 1);
  }

  hx = (xmax - xmin) / gnx;
  hy = (ymax - ymin) / gny;
  hz = (zmax - zmin) / gnz;
  global_index();
}

void init_vars(void)
{
  alpha = 1.0;
  beta = 0.0;
  init_KL();
}

void alloc_vars(void)
{
  int n;
  fieldgx = (int *)_mm_malloc(sizeof(int) * nx * ny * nz, 256);
  fieldgy = (int *)_mm_malloc(sizeof(int) * nx * ny * nz, 256);
  fieldgz = (int *)_mm_malloc(sizeof(int) * nx * ny * nz, 256);
  ac = (struct Allen_Cahn *)_mm_malloc(nac * sizeof(struct Allen_Cahn), 256);
  for (n = 0; n < nac; n++)
  {
    ac[n].fieldE = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldE1 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldE2 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldEt = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldEp = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldE1p = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldEs_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEr_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEs_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEr_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEs_top = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEr_top = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEs_bottom = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEr_bottom = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ac[n].fieldEr_front = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ac[n].fieldEr_back = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ac[n].fieldEe_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ac[n].fieldEe_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ac[n].fieldEe_top = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ac[n].fieldEe_bottom = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ac[n].fieldEe_front = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ac[n].fieldEe_back = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ac[n].fieldEu_left = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ac[n].fieldEu_right = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ac[n].fieldEu_top = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ac[n].fieldEu_bottom = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ac[n].fieldEu_front = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ac[n].fieldEu_back = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ac[n].fieldEmu_left = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ac[n].fieldEmu_right = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ac[n].fieldEmu_top = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ac[n].fieldEmu_bottom = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ac[n].fieldEmu_front = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ac[n].fieldEmu_back = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ac[n].phiE = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].phiE2 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].ireq_left_right_fieldE = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ac[n].ireq_top_bottom_fieldE = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ac[n].ireq_front_back_fieldE = (MPI_Request *)calloc(4, sizeof(MPI_Request));

    ac[n].felas = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].gradx = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].grady = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].gradz = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].f1 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].f2 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].f3 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
  }

  ch = (struct Cahn_Hilliard *)_mm_malloc(nch * sizeof(struct Cahn_Hilliard), 256);
  for (n = 0; n < nch; n++)
  {
    ch[n].f0 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].f1 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].f2 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].ft = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].M = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);

    ch[n].fieldCI = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCI1 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCI2 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCIt = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCIp = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCI1p = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCIs_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIr_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIs_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIr_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIs_top = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIr_top = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIs_bottom = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIr_bottom = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
    ch[n].fieldCIr_front = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ch[n].fieldCIr_back = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ch[n].fieldCIe_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ch[n].fieldCIe_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ch[n].fieldCIe_top = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ch[n].fieldCIe_bottom = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ch[n].fieldCIe_front = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ch[n].fieldCIe_back = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ch[n].fieldCIu_left = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ch[n].fieldCIu_right = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ch[n].fieldCIu_top = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ch[n].fieldCIu_bottom = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ch[n].fieldCIu_front = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ch[n].fieldCIu_back = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ch[n].fieldCImu_left = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ch[n].fieldCImu_right = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ch[n].fieldCImu_top = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ch[n].fieldCImu_bottom = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ch[n].fieldCImu_front = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ch[n].fieldCImu_back = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ch[n].phiCI = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].phiCI2 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].ireq_left_right_fieldCI = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ch[n].ireq_top_bottom_fieldCI = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ch[n].ireq_front_back_fieldCI = (MPI_Request *)calloc(4, sizeof(MPI_Request));

    ch[n].felas = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].felase_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ch[n].felase_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ch[n].felase_top = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ch[n].felase_bottom = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ch[n].felase_front = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ch[n].felase_back = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
  }
  MPX = (double *)_mm_malloc(sizeof(double) * nx * nx, 256);
  MPY = (double *)_mm_malloc(sizeof(double) * ny * ny, 256);
  MPZ = (double *)_mm_malloc(sizeof(double) * nz * nz, 256);
  MPXI = (double *)_mm_malloc(sizeof(double) * nx * nx, 256);
  MPYI = (double *)_mm_malloc(sizeof(double) * ny * ny, 256);
  MPZI = (double *)_mm_malloc(sizeof(double) * nz * nz, 256);
  DDX = (double *)_mm_malloc(sizeof(double) * nx, 256);
  DDY = (double *)_mm_malloc(sizeof(double) * ny, 256);
  DDZ = (double *)_mm_malloc(sizeof(double) * nz, 256);
  MPX_b = (double *)_mm_malloc(sizeof(double) * nx * nx * 4, 256);
  MPY_b = (double *)_mm_malloc(sizeof(double) * ny * ny * 4, 256);
  MPZ_b = (double *)_mm_malloc(sizeof(double) * nz * nz * 4, 256);
  MPXI_b = (double *)_mm_malloc(sizeof(double) * nx * nx * 4, 256);
  MPYI_b = (double *)_mm_malloc(sizeof(double) * ny * ny * 4, 256);
  MPZI_b = (double *)_mm_malloc(sizeof(double) * nz * nz * 4, 256);
  DDX_b = (double *)_mm_malloc(sizeof(double) * nx * 4, 256);
  DDY_b = (double *)_mm_malloc(sizeof(double) * ny * 4, 256);
  DDZ_b = (double *)_mm_malloc(sizeof(double) * nz * 4, 256);
  status = (MPI_Status *)calloc(4, sizeof(MPI_Status));
}

ac_init_phi(double LE, double KE, double *phi, double *phi2)
{
  int i, j, k, l;
  double Hijk, Gijk, tmp;
  for (j = iy1; j < iy4; j++)
  {
    for (i = ix1; i < ix4; i++)
    {
      for (k = iz1; k < iz4; k++)
      {
        l = j * nz * nx + i * nz + k;
        tmp = kkz * DDZ[k] + kky * DDY[j] + kkx * DDX[i];
        Gijk = tmp;
        Hijk = -LE * (tmp * epn2 - KE);
        if (fabs(Hijk) > 1.0e-8)
        {
          tmp = exp(-dt * Hijk);
          phi[l] = (1.0 - tmp) / Hijk;
          phi2[l] = (1.0 - phi[l] / dt) / Hijk;
        }
        else
        {
          phi[l] = dt;
          phi2[l] = dt / 2;
        }
      }
    }
  }
}

ch_init_phi(double LCI, double KCI, double *phi, double *phi2)
{
  int i, j, k, l;
  double Hijk, Gijk, tmp;
  for (j = iy1; j < iy4; j++)
  {
    for (i = ix1; i < ix4; i++)
    {
      for (k = iz1; k < iz4; k++)
      {
        l = j * nz * nx + i * nz + k;
        tmp = DDZ[k] + DDY[j] + DDX[i];
        Hijk = LCI * (tmp * tmp * epn2 - KCI * tmp);
        if (fabs(Hijk) > 1.0e-8)
        {
          tmp = exp(-dt * Hijk);
          phi[l] = (1.0 - tmp) / Hijk;
          phi2[l] = (1.0 - phi[l] / dt) / Hijk;
        }
        else
        {
          phi[l] = dt;
          phi2[l] = dt / 2;
        }
      }
    }
  }
}

void read_matrices(void)
{
  int n;
  char filename[1024];
  int i, j, k, l, id, mp_ofst, d_ofst;
  FILE *file;

  if (cart_id[1] == 0 && cart_id[2] == 0)
  {
    for (k = 0; k < 4; k++)
    {
      mp_ofst = k * nx * nx;
      d_ofst = k * nx;

      sprintf(filename, "d%d%d.dat", k, nx);
      file = fopen(filename, "r");
      for (l = 0; l < nx; l++)
      {
        fscanf(file, "%lf", DDX_b + d_ofst + l);
        DDX_b[d_ofst + l] = DDX_b[d_ofst + l];
      }
      fclose(file);

      sprintf(filename, "v%d%d.dat", k, nx);
      file = fopen(filename, "r");
      for (j = 0; j < nx; j++)
      {
        for (i = 0; i < nx; i++)
        {
          l = i * nx + j;
          fscanf(file, "%lf", MPX_b + mp_ofst + l);
        }
      }
      fclose(file);

      sprintf(filename, "vi%d%d.dat", k, nx);
      file = fopen(filename, "r");
      for (j = 0; j < nx; j++)
      {
        for (i = 0; i < nx; i++)
        {
          l = i * nx + j;
          fscanf(file, "%lf", MPXI_b + mp_ofst + l);
        }
      }
      fclose(file);

      mp_ofst = k * ny * ny;
      d_ofst = k * ny;

      sprintf(filename, "d%d%d.dat", k, ny);
      file = fopen(filename, "r");
      for (l = 0; l < ny; l++)
      {
        fscanf(file, "%lf", DDY_b + d_ofst + l);
      }
      fclose(file);

      sprintf(filename, "v%d%d.dat", k, ny);
      file = fopen(filename, "r");
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < ny; i++)
        {
          l = i * ny + j;
          fscanf(file, "%lf", MPY_b + mp_ofst + l);
        }
      }
      fclose(file);

      sprintf(filename, "vi%d%d.dat", k, ny);
      file = fopen(filename, "r");
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < ny; i++)
        {
          l = i * ny + j;
          fscanf(file, "%lf", MPYI_b + mp_ofst + l);
        }
      }
      fclose(file);

      mp_ofst = k * nz * nz;
      d_ofst = k * nz;

      sprintf(filename, "d%d%d.dat", k, nz);
      file = fopen(filename, "r");
      for (l = 0; l < nz; l++)
      {
        fscanf(file, "%lf", DDZ_b + d_ofst + l);
      }
      fclose(file);

      sprintf(filename, "v%d%d.dat", k, nz);
      file = fopen(filename, "r");
      for (j = 0; j < nz; j++)
      {
        for (i = 0; i < nz; i++)
        {
          l = i * nz + j;
          fscanf(file, "%lf", MPZ_b + mp_ofst + l);
        }
      }
      fclose(file);

      sprintf(filename, "vi%d%d.dat", k, nz);
      file = fopen(filename, "r");
      for (j = 0; j < nz; j++)
      {
        for (i = 0; i < nz; i++)
        {
          l = i * nz + j;
          fscanf(file, "%lf", MPZI_b + mp_ofst + l);
        }
      }
      fclose(file);
    }
  }
  MPI_Bcast(&DDX_b[0], nx * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&MPX_b[0], nx * nx * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&MPXI_b[0], nx * nx * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&DDY_b[0], ny * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&MPY_b[0], ny * ny * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&MPYI_b[0], ny * ny * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&DDZ_b[0], nz * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&MPZ_b[0], nz * nz * 4, MPI_DOUBLE, 0, YZ_COMM);
  MPI_Bcast(&MPZI_b[0], nz * nz * 4, MPI_DOUBLE, 0, YZ_COMM);
  if (left < 0)
  {
    k = 1;
  }
  else
  {
    k = 2;
  }

  if (right < 0)
  {
    l = 1;
  }
  else
  {
    l = 2;
  }

  id = (k - 1) + (l - 1) * 2;
  mp_ofst = id * nx * nx;
  d_ofst = id * nx;

  for (l = 0; l < nx; l++)
  {
    DDX[l] = DDX_b[d_ofst + l] / hx / hx;
  }
  for (l = 0; l < nx * nx; l++)
  {
    MPX[l] = MPX_b[mp_ofst + l];
    MPXI[l] = MPXI_b[mp_ofst + l];
  }

  if (top < 0)
  {
    k = 1;
  }
  else
  {
    k = 2;
  }

  if (bottom < 0)
  {
    l = 1;
  }
  else
  {
    l = 2;
  }

  id = (k - 1) + (l - 1) * 2;
  mp_ofst = id * ny * ny;
  d_ofst = id * ny;

  for (l = 0; l < ny; l++)
  {
    DDY[l] = DDY_b[d_ofst + l] / hy / hy;
  }
  for (l = 0; l < ny * ny; l++)
  {
    MPY[l] = MPY_b[mp_ofst + l];
    MPYI[l] = MPYI_b[mp_ofst + l];
  }

  if (front < 0)
  {
    k = 1;
  }
  else
  {
    k = 2;
  }

  if (back < 0)
  {
    l = 1;
  }
  else
  {
    l = 2;
  }

  id = (k - 1) + (l - 1) * 2;
  mp_ofst = id * nz * nz;
  d_ofst = id * nz;

  for (l = 0; l < nz; l++)
  {
    DDZ[l] = DDZ_b[d_ofst + l] / hz / hz;
  }
  for (l = 0; l < nz * nz; l++)
  {
    MPZ[l] = MPZ_b[mp_ofst + l];
    MPZI[l] = MPZI_b[mp_ofst + l];
  }
  for (n = 0; n < nac; n++)
  {
    ac_init_phi(ac[n].LE, ac[n].KE, ac[n].phiE, ac[n].phiE2);
  }
  for (n = 0; n < nch; n++)
  {
    ch_init_phi(ch[n].LCI, ch[n].KCI, ch[n].phiCI, ch[n].phiCI2);
  }
}

void dealloc_vars(void)
{
  int l;
  int n;
  _mm_free(fieldgx);
  _mm_free(fieldgy);
  _mm_free(fieldgz);
  for (n = 0; n < nac; n++)
  {
    _mm_free(ac[n].fieldE);
    _mm_free(ac[n].fieldE1);
    _mm_free(ac[n].fieldE2);
    _mm_free(ac[n].fieldEt);
    _mm_free(ac[n].fieldEp);
    _mm_free(ac[n].fieldE1p);
    _mm_free(ac[n].fieldEs_left);
    _mm_free(ac[n].fieldEr_left);
    _mm_free(ac[n].fieldEs_right);
    _mm_free(ac[n].fieldEr_right);
    _mm_free(ac[n].fieldEs_top);
    _mm_free(ac[n].fieldEr_top);
    _mm_free(ac[n].fieldEs_bottom);
    _mm_free(ac[n].fieldEr_bottom);
    _mm_free(ac[n].fieldEr_front);
    _mm_free(ac[n].fieldEr_back);
    _mm_free(ac[n].fieldEe_left);
    _mm_free(ac[n].fieldEe_right);
    _mm_free(ac[n].fieldEe_top);
    _mm_free(ac[n].fieldEe_bottom);
    _mm_free(ac[n].fieldEe_front);
    _mm_free(ac[n].fieldEe_back);
    _mm_free(ac[n].fieldEu_left);
    _mm_free(ac[n].fieldEu_right);
    _mm_free(ac[n].fieldEu_top);
    _mm_free(ac[n].fieldEu_bottom);
    _mm_free(ac[n].fieldEu_front);
    _mm_free(ac[n].fieldEu_back);
    _mm_free(ac[n].fieldEmu_left);
    _mm_free(ac[n].fieldEmu_right);
    _mm_free(ac[n].fieldEmu_top);
    _mm_free(ac[n].fieldEmu_bottom);
    _mm_free(ac[n].fieldEmu_front);
    _mm_free(ac[n].fieldEmu_back);
    _mm_free(ac[n].phiE);
    _mm_free(ac[n].phiE2);
    for (l = 0; l < 4; l++)
    {
      MPI_Request_free(&ac[n].ireq_left_right_fieldE[l]);
      MPI_Request_free(&ac[n].ireq_top_bottom_fieldE[l]);
      MPI_Request_free(&ac[n].ireq_front_back_fieldE[l]);
    }
    free(ac[n].ireq_left_right_fieldE);
    free(ac[n].ireq_top_bottom_fieldE);
    free(ac[n].ireq_front_back_fieldE);
    _mm_free(ac[n].felas);
    _mm_free(ac[n].gradx);
    _mm_free(ac[n].grady);
    _mm_free(ac[n].gradz);
    _mm_free(ac[n].f1);
    _mm_free(ac[n].f2);
    _mm_free(ac[n].f3);
  }
  _mm_free(ac);

  for (n = 0; n < nch; n++)
  {
    _mm_free(ch[n].fieldCI);
    _mm_free(ch[n].fieldCI1);
    _mm_free(ch[n].fieldCI2);
    _mm_free(ch[n].fieldCIt);
    _mm_free(ch[n].fieldCIp);
    _mm_free(ch[n].fieldCI1p);
    _mm_free(ch[n].fieldCIs_left);
    _mm_free(ch[n].fieldCIr_left);
    _mm_free(ch[n].fieldCIs_right);
    _mm_free(ch[n].fieldCIr_right);
    _mm_free(ch[n].fieldCIs_top);
    _mm_free(ch[n].fieldCIr_top);
    _mm_free(ch[n].fieldCIs_bottom);
    _mm_free(ch[n].fieldCIr_bottom);
    _mm_free(ch[n].fieldCIr_front);
    _mm_free(ch[n].fieldCIr_back);
    _mm_free(ch[n].fieldCIe_left);
    _mm_free(ch[n].fieldCIe_right);
    _mm_free(ch[n].fieldCIe_top);
    _mm_free(ch[n].fieldCIe_bottom);
    _mm_free(ch[n].fieldCIe_front);
    _mm_free(ch[n].fieldCIe_back);
    _mm_free(ch[n].fieldCIu_left);
    _mm_free(ch[n].fieldCIu_right);
    _mm_free(ch[n].fieldCIu_top);
    _mm_free(ch[n].fieldCIu_bottom);
    _mm_free(ch[n].fieldCIu_front);
    _mm_free(ch[n].fieldCIu_back);
    _mm_free(ch[n].fieldCImu_left);
    _mm_free(ch[n].fieldCImu_right);
    _mm_free(ch[n].fieldCImu_top);
    _mm_free(ch[n].fieldCImu_bottom);
    _mm_free(ch[n].fieldCImu_front);
    _mm_free(ch[n].fieldCImu_back);
    _mm_free(ch[n].phiCI);
    _mm_free(ch[n].phiCI2);
    for (l = 0; l < 4; l++)
    {
      MPI_Request_free(&ch[n].ireq_left_right_fieldCI[l]);
      MPI_Request_free(&ch[n].ireq_top_bottom_fieldCI[l]);
      MPI_Request_free(&ch[n].ireq_front_back_fieldCI[l]);
    }
    free(ch[n].ireq_left_right_fieldCI);
    free(ch[n].ireq_top_bottom_fieldCI);
    free(ch[n].ireq_front_back_fieldCI);
    _mm_free(ch[n].felas);
    _mm_free(ch[n].felase_left);
    _mm_free(ch[n].felase_right);
    _mm_free(ch[n].felase_top);
    _mm_free(ch[n].felase_bottom);
    _mm_free(ch[n].felase_front);
    _mm_free(ch[n].felase_back);
  }
  _mm_free(ch);
  _mm_free(MPX);
  _mm_free(MPY);
  _mm_free(MPZ);
  _mm_free(MPXI);
  _mm_free(MPYI);
  _mm_free(MPZI);
  _mm_free(DDX);
  _mm_free(DDY);
  _mm_free(DDZ);
  _mm_free(MPX_b);
  _mm_free(MPY_b);
  _mm_free(MPZ_b);
  _mm_free(MPXI_b);
  _mm_free(MPYI_b);
  _mm_free(MPZI_b);
  _mm_free(DDX_b);
  _mm_free(DDY_b);
  _mm_free(DDZ_b);
  free(status);
}
