#include "ScLETD.h"
#include "fft3d.h"
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif
#define DIM 3
int NX, NY, NZ;
// [11, 22, 33, 23, 13, 12]
#define NELASTIC 4

double epsilon2d[NELASTIC][DIM * DIM];
double sigma2d[NELASTIC][DIM * DIM];
double C4d[DIM * DIM * DIM * DIM];
double *Bn[NELASTIC][NELASTIC];
double *f;
double *yy, *elas;
double *Bn_starre, *Bn_starim, *Bn_stariftre;
double *tmpy_re[NELASTIC], *tmpy_fftre, *tmpy_fftim;
double *s_left, *s_right, *s_top, *s_bottom, *s_front, *s_back;
double *r_left, *r_right, *r_top, *r_bottom, *r_front, *r_back;
double *e_left, *e_right, *e_top, *e_bottom, *e_front, *e_back;
MPI_Request *ireq_left_right, *ireq_top_bottom, *ireq_front_back;

static void
write_small(double *f)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  // z section
  sprintf(fname, "./Bn/Bnz%d%d%d", cart_id[0], cart_id[1], cart_id[2]);
  //if (cart_id[2] == (procs[2]/2)) {
  if (cart_id[2] == 0)
  {
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
      printf("fopen error %s!\n", strerror(errno));
      exit(1);
    }
    k = 0;
    for (j = 0; j < NY; j++)
    {
      for (i = 0; i < NX; i++)
      {
        fprintf(fp, "%+1.15lf ", f[k * NX * NY + j * NX + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  // x section
  sprintf(fname, "./Bn/Bnx%d%d%d", cart_id[1], cart_id[2], cart_id[0]);
  if (cart_id[0] == (procs[0] / 2))
  {
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
      printf("fopen error %s!\n", strerror(errno));
      exit(1);
    }
    i = 0;
    for (k = 0; k < NZ; k++)
    {
      for (j = 0; j < NY; j++)
      {
        fprintf(fp, "%+1.15lf ", f[k * NX * NY + j * NX + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  sprintf(fname, "./Bn/Bny%d%d%d", cart_id[0], cart_id[2], cart_id[1]);
  if (cart_id[1] == (procs[1] / 2))
  {
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
      printf("fopen error %s!\n", strerror(errno));
      exit(1);
    }
    j = 0;
    for (k = 0; k < NZ; k++)
    {
      for (i = 0; i < NX; i++)
      {
        fprintf(fp, "%+1.15lf ", f[k * NX * NY + j * NX + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

static void
write_large(char *str, double *f)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  // z section
  sprintf(fname, "./output/%s_z%d%d%d", str, cart_id[0], cart_id[1], cart_id[2]);
  //if (cart_id[2] == (procs[2]/2)) {
  if (cart_id[2] == 0)
  {
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
      printf("fopen error %s!\n", strerror(errno));
      exit(1);
    }
    k = 6;
    for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
      {
        fprintf(fp, "%+1.15lf ", f[k * nx * ny + j * nx + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  // x section
  sprintf(fname, "./output/%s_x%d%d%d", str, cart_id[1], cart_id[2], cart_id[0]);
  if (cart_id[0] == (procs[0] / 2))
  {
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
      printf("fopen error %s!\n", strerror(errno));
      exit(1);
    }
    i = 6;
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        fprintf(fp, "%+1.15lf ", f[k * nx * ny + j * nx + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  sprintf(fname, "./output/%s_y%d%d%d", str, cart_id[0], cart_id[2], cart_id[1]);
  if (cart_id[1] == (procs[1] / 2))
  {
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
      printf("fopen error %s!\n", strerror(errno));
      exit(1);
    }
    j = 6;
    for (k = 0; k < nz; k++)
    {
      for (i = 0; i < nx; i++)
      {
        fprintf(fp, "%+1.15lf ", f[k * nx * ny + j * nx + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

static void
check_max(double *v)
{
  int i, j, k;
  double tmp, maxtmp, mintmp, mmax, mmin;
  maxtmp = -1.0e30;
  mintmp = 1.0e30;
  for (k = 0; k < NZ; k++)
  {
    for (j = 0; j < NY; j++)
    {
      for (i = 0; i < NX; i++)
      {
        tmp = v[k * NX * NY + j * NX + i];
        if (tmp > maxtmp)
        {
          maxtmp = tmp;
        }
        if (tmp < mintmp)
        {
          mintmp = tmp;
        }
      }
    }
  }
  MPI_Reduce(&maxtmp, &mmax, 1, MPI_DOUBLE, MPI_MAX, prank, MPI_COMM_WORLD);
  MPI_Reduce(&mintmp, &mmin, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);

  if (myrank == prank)
  {
    printf("max:\t\t%+1.15lf,\t\t%+1.15lf\n", mmax, mmin);
  }
}

// local variables
static void
elastic_malloc()
{
  int p, q;
  for (p = 0; p < NELASTIC; p++)
  {
    for (q = 0; q < NELASTIC; q++)
    {
      Bn[p][q] = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
    }
    tmpy_re[p] = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  }
  Bn_starre = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  Bn_starim = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  Bn_stariftre = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  tmpy_fftre = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  tmpy_fftim = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  yy = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);
  elas = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ, 256);

  f = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
  s_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
  s_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
  s_top = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
  s_bottom = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
  s_front = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
  s_back = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
  r_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
  r_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2), 256);
  r_top = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
  r_bottom = (double *)_mm_malloc(sizeof(double) * nx * (nghost + 2) * (nz + (nghost + 2) * 2), 256);
  r_front = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
  r_back = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
  e_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
  e_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
  e_top = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
  e_bottom = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
  e_front = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
  e_back = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
  ireq_left_right = (MPI_Request *)calloc(4, sizeof(MPI_Request));
  ireq_top_bottom = (MPI_Request *)calloc(4, sizeof(MPI_Request));
  ireq_front_back = (MPI_Request *)calloc(4, sizeof(MPI_Request));
}

static void
elastic_init_transfer()
{
  MPI_Send_init(&s_left[0], 1, left_right, left, 9, MPI_COMM_WORLD, &ireq_left_right[0]);
  MPI_Recv_init(&r_right[0], 1, left_right, right, 9, MPI_COMM_WORLD, &ireq_left_right[1]);
  MPI_Send_init(&s_right[0], 1, left_right, right, 9, MPI_COMM_WORLD, &ireq_left_right[2]);
  MPI_Recv_init(&r_left[0], 1, left_right, left, 9, MPI_COMM_WORLD, &ireq_left_right[3]);
  MPI_Send_init(&s_top[0], 1, top_bottom, top, 9, MPI_COMM_WORLD, &ireq_top_bottom[0]);
  MPI_Recv_init(&r_bottom[0], 1, top_bottom, bottom, 9, MPI_COMM_WORLD, &ireq_top_bottom[1]);
  MPI_Send_init(&s_bottom[0], 1, top_bottom, bottom, 9, MPI_COMM_WORLD, &ireq_top_bottom[2]);
  MPI_Recv_init(&r_top[0], 1, top_bottom, top, 9, MPI_COMM_WORLD, &ireq_top_bottom[3]);
  MPI_Send_init(&f[nx * ny * nghost], 1, front_back, front, 9, MPI_COMM_WORLD, &ireq_front_back[0]);
  MPI_Recv_init(&r_back[0], 1, front_back, back, 9, MPI_COMM_WORLD, &ireq_front_back[1]);
  MPI_Send_init(&f[nx * ny * (nz - nghost - (nghost + 2))], 1, front_back, back, 9, MPI_COMM_WORLD, &ireq_front_back[2]);
  MPI_Recv_init(&r_front[0], 1, front_back, front, 9, MPI_COMM_WORLD, &ireq_front_back[3]);
}

void elastic_finish()
{
  int p, q, l;
  fft_finish();
  for (p = 0; p < NELASTIC; p++)
  {
    for (q = 0; q < NELASTIC; q++)
    {
      _mm_free(Bn[p][q]);
    }
    _mm_free(tmpy_re[p]);
  }
  _mm_free(Bn_starre);
  _mm_free(Bn_starim);
  _mm_free(Bn_stariftre);
  _mm_free(tmpy_fftre);
  _mm_free(tmpy_fftim);
  _mm_free(yy);
  _mm_free(elas);

  /* free elastic variable*/
  _mm_free(f);
  _mm_free(s_left);
  _mm_free(s_right);
  _mm_free(s_top);
  _mm_free(s_bottom);
  _mm_free(s_front);
  _mm_free(s_back);
  _mm_free(r_left);
  _mm_free(r_right);
  _mm_free(r_top);
  _mm_free(r_bottom);
  _mm_free(r_front);
  _mm_free(r_back);
  _mm_free(e_left);
  _mm_free(e_right);
  _mm_free(e_top);
  _mm_free(e_bottom);
  _mm_free(e_front);
  _mm_free(e_back);
  for (l = 0; l < 4; l++)
  {
    MPI_Request_free(&ireq_left_right[l]);
    MPI_Request_free(&ireq_top_bottom[l]);
    MPI_Request_free(&ireq_front_back[l]);
  }
  free(ireq_left_right);
  free(ireq_top_bottom);
  free(ireq_front_back);
}

static void
elastic_index_3to6(const int i, const int j, int *m)
{
  if (i == 0 && j == 0)
  {
    *m = 0;
  }
  else if (i == 1 && j == 1)
  {
    *m = 1;
  }
  else if (i == 2 && j == 2)
  {
    *m = 2;
  }
  else if ((i == 1 && j == 2) || (i == 2 && j == 1))
  {
    *m = 3;
  }
  else if ((i == 0 && j == 2) || (i == 2 && j == 0))
  {
    *m = 4;
  }
  else if ((i == 0 && j == 1) || (i == 1 && j == 0))
  {
    *m = 5;
  }
}

void elastic_input()
{

  int i, j, k, l;
  int m, n;
  FILE *fp;
  char *aline = NULL;
  size_t len = 0;
  double C2d[6 * 6];
  NX = nx - 2 * nghost;
  NY = ny - 2 * nghost;
  NZ = nz - 2 * nghost;

  fp = fopen("elastic_input.txt", "r");
  if (fp == NULL)
  {
    printf("open elastic_input.txt Fail!");
    exit(EXIT_FAILURE);
  }
  if (myrank == 0)
  {
    printf("Reading input.txt file:\n");
  }
  //getline(&aline, &len, fp);
  //printf("%s\n", aline);

  fscanf(fp, "%lf\n", &ElasticScale);
  /*
  C00 C01 C02 C03 C04 C05 
      C11 C12 C12 C14 C15
          C22 C23 C24 C25
              C33 C34 C35
                  C44 C45
                      C55
  */
  for (i = 0; i < 2 * DIM; i++)
  {
    for (j = 0; j < 2 * DIM; j++)
    {
      fscanf(fp, "%lf\n", &C2d[i * 2 * DIM + j]);
    }
  }

  // 00: 0, 11: 1, 22: 2, 12/21: 3, 02/20: 4, 01/10: 5
  for (i = 0; i < DIM; i++)
  {
    for (j = 0; j < DIM; j++)
    {
      for (k = 0; k < DIM; k++)
      {
        for (l = 0; l < DIM; l++)
        {
          elastic_index_3to6(i, j, &m);
          elastic_index_3to6(k, l, &n);
          C4d[i * DIM * DIM * DIM + j * DIM * DIM + k * DIM + l] = C2d[m * 2 * DIM + n];
        }
      }
    }
  }
  int p;
  for (p = 0; p < NELASTIC; p++)
  {
    // 00: 0, 11: 1, 22: 2, 12/21: 3, 02/20: 4, 01/10: 5
    for (i = 0; i < DIM; i++)
    {
      for (j = 0; j < DIM; j++)
      {
        fscanf(fp, "%lf", &(epsilon2d[p][i * DIM + j]));
      }
    }
    // sigma0[ij] = C[ijkl] * epsilon0[kl]
    for (i = 0; i < DIM; i++)
    {
      for (j = 0; j < DIM; j++)
      {
        sigma2d[p][i * DIM + j] = 0;
        for (k = 0; k < DIM; k++)
        {
          for (l = 0; l < DIM; l++)
          {
            sigma2d[p][i * DIM + j] += C4d[i * DIM * DIM * DIM + j * DIM * DIM + k * DIM + l] * epsilon2d[p][k * DIM + l];
          }
        }
      }
    }
  }
  fclose(fp);
  if (aline)
    free(aline);
}

// omega = 1/normal(omega_inver) * (-1)^(i+j) * M(omega_inver)
static void
inverse_3x3(double *omega, double *omega_inver)
{
  double normal = 0.0;
  normal += omega_inver[0 * DIM + 0] * (omega_inver[1 * DIM + 1] * omega_inver[2 * DIM + 2] - omega_inver[1 * DIM + 2] * omega_inver[2 * DIM + 1]);
  normal -= omega_inver[0 * DIM + 1] * (omega_inver[1 * DIM + 0] * omega_inver[2 * DIM + 2] - omega_inver[1 * DIM + 2] * omega_inver[2 * DIM + 0]);
  normal += omega_inver[0 * DIM + 2] * (omega_inver[1 * DIM + 0] * omega_inver[2 * DIM + 1] - omega_inver[1 * DIM + 1] * omega_inver[2 * DIM + 0]);
  omega[0 * DIM + 0] = (+1) * (omega_inver[1 * DIM + 1] * omega_inver[2 * DIM + 2] - omega_inver[1 * DIM + 2] * omega_inver[2 * DIM + 1]) / normal;
  omega[0 * DIM + 1] = (-1) * (omega_inver[1 * DIM + 0] * omega_inver[2 * DIM + 2] - omega_inver[2 * DIM + 0] * omega_inver[1 * DIM + 2]) / normal;
  omega[0 * DIM + 2] = (+1) * (omega_inver[1 * DIM + 0] * omega_inver[2 * DIM + 1] - omega_inver[2 * DIM + 0] * omega_inver[1 * DIM + 1]) / normal;
  omega[1 * DIM + 0] = (-1) * (omega_inver[0 * DIM + 1] * omega_inver[2 * DIM + 2] - omega_inver[2 * DIM + 1] * omega_inver[0 * DIM + 2]) / normal;
  omega[1 * DIM + 1] = (+1) * (omega_inver[0 * DIM + 0] * omega_inver[2 * DIM + 2] - omega_inver[0 * DIM + 2] * omega_inver[2 * DIM + 0]) / normal;
  omega[1 * DIM + 2] = (-1) * (omega_inver[0 * DIM + 0] * omega_inver[2 * DIM + 1] - omega_inver[2 * DIM + 0] * omega_inver[0 * DIM + 1]) / normal;
  omega[2 * DIM + 0] = (+1) * (omega_inver[0 * DIM + 1] * omega_inver[1 * DIM + 2] - omega_inver[1 * DIM + 1] * omega_inver[0 * DIM + 2]) / normal;
  omega[2 * DIM + 1] = (-1) * (omega_inver[0 * DIM + 0] * omega_inver[1 * DIM + 2] - omega_inver[1 * DIM + 0] * omega_inver[0 * DIM + 2]) / normal;
  omega[2 * DIM + 2] = (+1) * (omega_inver[0 * DIM + 0] * omega_inver[1 * DIM + 1] - omega_inver[0 * DIM + 1] * omega_inver[1 * DIM + 0]) / normal;
}

static void
elastic_calculate_BN()
{
  int p, q;
  int x, y, z;
  int i, j, k, l;
  double n11[NX], n22[NY], n33[NZ];
  double normal;

  // calculate n
  int gnnx = procs[0] * NX; // 256
  int gnny = procs[1] * NY; // 256
  int gnnz = procs[2] * NZ; // 256
  int cntx = gnnx / 2;      // 128
  int cnty = gnny / 2;      // 128
  int cntz = gnnz / 2;      // 128
  for (x = 0; x < NX; x++)  // 0-64
  {
    int gx = cart_id[0] * NX + x; // 0-255
    if (gx < cntx)                // 0-127
    {
      n11[x] = 1.0 * gx; // 0-127
    }
    else
    {                              // 128-255
      n11[x] = -1.0 * (gnnx - gx); // (-128)-(-1)
    }
  }

  for (y = 0; y < NY; y++) // 0-64
  {
    int gy = cart_id[1] * NY + y; // 0-255
    if (gy < cnty)                // 0-127
    {
      n22[y] = 1.0 * gy; // 0-127
    }
    else
    {                              // 128-255
      n22[y] = -1.0 * (gnny - gy); // (-128)-(-1)
    }
  }
  for (z = 0; z < NZ; z++) // 0-64
  {
    int gz = cart_id[2] * NZ + z; // 0-255
    if (gz < cntz)                // 0-127
    {
      n33[z] = 1.0 * gz; // 0-127
    }
    else
    {                              // 128-255
      n33[z] = -1.0 * (gnnz - gz); // (-128)-(-1)
    }
  }
  /* 
    Bn(p,q) = C[ijkl] * epsilon0(p)[ij] * epsilon0(q)[kl] 
     - n[i] * sigma0(p)[ij] * omega[jk] * sigma0(q)[kl] * n[l]
  */
  for (p = 0; p < NELASTIC; p++)
  {
    for (q = 0; q < NELASTIC; q++)
    {
      for (z = 0; z < NZ; z++)
      {
        for (y = 0; y < NY; y++)
        {
          for (x = 0; x < NX; x++)
          {
            double n123[DIM];
            double omega_inver[DIM * DIM];
            double omega[DIM * DIM];
            double tmp0, tmp1;

            int gx = cart_id[0] * NX + x;
            int gy = cart_id[1] * NY + y;
            int gz = cart_id[2] * NZ + z;
            if (gx == 0 && gy == 0 && gz == 0)
            {
              n123[0] = 0;
              n123[1] = 1;
              n123[2] = 0;
            }
            else
            {
              n123[0] = n11[x] / pow((pow(n11[x], 2) + pow(n22[y], 2) + pow(n33[z], 2)), 0.5);
              n123[1] = n22[y] / pow((pow(n11[x], 2) + pow(n22[y], 2) + pow(n33[z], 2)), 0.5);
              n123[2] = n33[z] / pow((pow(n11[x], 2) + pow(n22[y], 2) + pow(n33[z], 2)), 0.5);
            }

            // omega_inver[ik] = C[ijkl] * n[j] * n[l]
            for (i = 0; i < DIM; ++i)
            {
              for (k = 0; k < DIM; ++k)
              {
                omega_inver[i * DIM + k] = 0.0;
                for (j = 0; j < DIM; ++j)
                {
                  for (l = 0; l < DIM; ++l)
                  {
                    omega_inver[i * DIM + k] += C4d[i * DIM * DIM * DIM + j * DIM * DIM + k * DIM + l] * n123[j] * n123[l];
                  }
                }
              }
            }

            // omega = inverse(omega_inver)
            inverse_3x3(omega, omega_inver);

            tmp0 = 0;
            // tmp0 = C[ijkl] * epsilon[ij] * epsilon[kl]
            for (i = 0; i < DIM; i++)
            {
              for (j = 0; j < DIM; j++)
              {
                for (k = 0; k < DIM; k++)
                {
                  for (l = 0; l < DIM; l++)
                  {
                    tmp0 += C4d[i * DIM * DIM * DIM + j * DIM * DIM + k * DIM + l] * epsilon2d[p][i * DIM + j] * epsilon2d[q][k * DIM + l];
                  }
                }
              }
            }

            // tmp1 = n[i] * sigma0[ij] * omega[jk] * sigma0[kl] *n[l]
            tmp1 = 0;
            for (i = 0; i < DIM; i++)
            {
              for (j = 0; j < DIM; j++)
              {
                for (k = 0; k < DIM; k++)
                {
                  for (l = 0; l < DIM; l++)
                  {
                    tmp1 += n123[i] * sigma2d[p][i * DIM + j] * omega[j * DIM + k] * sigma2d[q][k * DIM + l] * n123[l];
                  }
                }
              }
            }
            // Bn = tmp0 - tmp1
            Bn[p][q][z * NY * NX + y * NX + x] = tmp0 - tmp1;
            //Bn[z * NY * NX + y * NX + x] = 75 - tmp1;
          }
        }
      }
    }
  }
}

void elastic_init()
{
  fft_setup();
  elastic_malloc();
  elastic_init_transfer();
  elastic_calculate_BN();
}

// elastic = dtheta(eta(p)) * IFFT(Bn(p,q) * FFT(theta(eta(q))))
static void
elastic_calculate_ElasDri(int p)
{
  int x, y, z;
  int q;
  for (x = 0; x < NX; x++)
  {
    for (y = 0; y < NY; y++)
    {
      for (z = 0; z < NZ; z++)
      {
        elas[z * NY * NX + y * NX + x] = 0;
      }
    }
  }
  for (q = 0; q < NELASTIC; q++)
  {
    // tmpy_fft = FFT(theta(c))
    fft_forward(tmpy_re[q], tmpy_fftre, tmpy_fftim);

    // Bn_star = Bn * FFT(tmpy_fft)
    for (x = 0; x < NX; x++)
    {
      for (y = 0; y < NY; y++)
      {
        for (z = 0; z < NZ; z++)
        {
          Bn_starre[z * NY * NX + y * NX + x] = tmpy_fftre[z * NY * NX + y * NX + x] * Bn[p][q][z * NY * NX + y * NX + x];
          Bn_starim[z * NY * NX + y * NX + x] = tmpy_fftim[z * NY * NX + y * NX + x] * Bn[p][q][z * NY * NX + y * NX + x];
        }
      }
    }

    // Bn_starift = IFFT(Bn_star)
    fft_backward(Bn_starre, Bn_starim, Bn_stariftre);
    // elas = 2 * c * Bn_stariftre
    for (x = 0; x < NX; x++)
    {
      for (y = 0; y < NY; y++)
      {
        for (z = 0; z < NZ; z++)
        {
          elas[z * NY * NX + y * NX + x] += Bn_stariftre[z * NY * NX + y * NX + x];
        }
      }
    }
  }

  // transfer and get result in f
  for (x = 0; x < NX; x++)
  {
    for (y = 0; y < NY; y++)
    {
      for (z = 0; z < NZ; z++)
      {
        f[(z + nghost) * ny * nx + (y + nghost) * nx + x + nghost] =
            elas[z * NY * NX + y * NX + x];
      }
    }
  }
  MPI_Startall(4, ireq_front_back);
  MPI_Waitall(4, ireq_front_back, status);

  top_bottom_pack(f, s_top, s_bottom, r_front, r_back);
  MPI_Startall(4, ireq_top_bottom);
  MPI_Waitall(4, ireq_top_bottom, status);

  left_right_pack(f, s_left, s_right, r_top, r_bottom, r_front, r_back);
  MPI_Startall(4, ireq_left_right);
  MPI_Waitall(4, ireq_left_right, status);
  unpack(f, r_left, r_right, r_top, r_bottom, r_front, r_back);
}

void elastic_calculate()
{
  int p;
  int x, y, z;
  // copy all field[p]
  if (ELASTIC == AC_FUNCTION) {
    for (p = 0; p < nac; p++)
    {
      for (x = 0; x < NX; x++)
      {
        for (y = 0; y < NY; y++)
        {
          for (z = 0; z < NZ; z++)
          {
            // theta(eta) = eta, dtheta(eta) = 1
            tmpy_re[p][z * NY * NX + y * NX + x] = ac[p].fieldE[(z + nghost) * ny * nx + (y + nghost) * nx + x + nghost];
          }
        }
      }
    }
    for (p = 0; p < nac; p++)
    {
      // calculate elas of field[p]
      elastic_calculate_ElasDri(p);
      // copy elastic back to function
      for (x = 0; x < nx; x++)
      {
        for (y = 0; y < ny; y++)
        {
          for (z = 0; z < nz; z++)
          {
            ac[p].felas[z * ny * nx + y * nx + x] = f[z * ny * nx + y * nx + x];
          }
        }
      }
    }
  }
  if (ELASTIC == CH_FUNCTION) {
    for (p = 0; p < nch; p++)
    {
      for (x = 0; x < NX; x++)
      {
        for (y = 0; y < NY; y++)
        {
          for (z = 0; z < NZ; z++)
          {
            // theta(eta) = eta, dtheta(eta) = 1
            tmpy_re[p][z * NY * NX + y * NX + x] = 
              ch[p].fieldCI[(z + nghost) * ny * nx + (y + nghost) * nx + x + nghost];
          }
        }
      }
    }
    for (p = 0; p < nch; p++)
    {
      // calculate elas of field[p]
      elastic_calculate_ElasDri(p);
      // copy elastic back to function
      for (x = 0; x < nx; x++)
      {
        for (y = 0; y < ny; y++)
        {
          for (z = 0; z < nz; z++)
          {
            ch[p].felas[z * ny * nx + y * nx + x] = f[z * ny * nx + y * nx + x];
          }
        }
      }
      enlarge(ch[p].felase_left, ch[p].felase_right, ch[p].felase_top, 
          ch[p].felase_bottom, ch[p].felase_front, ch[p].felase_back, 
          r_left, r_right, r_top, r_bottom, r_front, r_back);
    }
  }
}

