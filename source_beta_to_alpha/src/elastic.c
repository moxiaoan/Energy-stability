#include "ScLETD.h"
#include "fft3d.h"
#include <string.h>
#include <errno.h>
#include <hipfft.h>
#include <math.h>
#include <stdio.h>
//#include "hip_complex.h"
#include <stdlib.h>
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif
#include <time.h>
void
write_small(double *f)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  // z section
  sprintf(fname, "./Bn/Bnz%d%d%d", cart_id[0], cart_id[1], cart_id[2]);
  if (cart_id[2] == (procs[2]/2)) {
  //if (cart_id[2] == 0)
//  {
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
void
//write_n(double *n, int a, int b)
write_n(double *n)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  sprintf(fname, "./data/omega_%d%d%d", cart_id[0], cart_id[1], cart_id[2]);
  //sprintf(fname, "./data/n11_%d%d%d", cart_id[0], cart_id[1], cart_id[2]);
  fp = fopen (fname, "w");
  for (k = 0; k < DIM; k++) {
  	for (j = 0; j < DIM; j++) {
    		fprintf (fp, "%+1.15lf\n", n[k*DIM+j]);
  	}
    		//fprintf (fp, "%+1.15lf\n", n[k]);
  }
  fclose (fp);
}
void
write_Bn_re(double *f, int p, int q)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  sprintf(fname, "./Bn/Bn_%d%d_%d%d%d", p, q, cart_id[0], cart_id[1], cart_id[2]);
  fp = fopen (fname, "w");
  for (k = 0; k < 16; k++) {
    for (j = 0; j < 16; j++) {
      for (i = 0; i < 16; i++) {
        fprintf (fp, "%+1.15lf\n", f[k * 16 * 16 + j * 16 + i]);
      }
    }
  }
  fclose (fp);
}
void
write_Bn_im(double *f, int p)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  sprintf(fname, "./Bn/Bn_im_%d_%d%d%d", p, cart_id[0], cart_id[1], cart_id[2]);
  fp = fopen (fname, "w");
  for (k = 0; k < NZ; k++) {
    for (j = 0; j < NY; j++) {
      for (i = 0; i < NX; i++) {
        fprintf (fp, "%+1.15lf\n", f[k * NX * NY + j * NX + i]);
      }
    }
  }
  fclose (fp);
}

void
cal_occupy_im_re(double *im, double *re)
{
  int i, j, k;
  double Re, Im;
  double maxtmp_re, maxtmp_im, mmax_re, mmax_im;
  maxtmp_re = -1.0e30;
  maxtmp_im = -1.0e30;
  for (k = 0; k < NZ; k++) {
    for (j = 0; j < NY; j++) {
      for (i = 0; i < NX; i++) {
        Im = im[k * NX * NY + j * NX + i];
        Re = re[k * NX * NY + j * NX + i];
        if (fabs(Re) > maxtmp_re)
        {
          maxtmp_re = fabs(Re);
        }
        if (fabs(Im) > maxtmp_im)
        {
          maxtmp_im = fabs(Im);
        }
      }
    }
  }
  MPI_Reduce(&maxtmp_re, &mmax_re, 1, MPI_DOUBLE, MPI_MAX, prank, MPI_COMM_WORLD);
  MPI_Reduce(&maxtmp_im, &mmax_im, 1, MPI_DOUBLE, MPI_MAX, prank, MPI_COMM_WORLD);

  if (myrank == prank)
  {
    printf("Re:\t\t%+1.15lf, Im:\t\t%+1.15lf, Ratio=Im/Re:\t\t%+1.15lf\n", mmax_re, mmax_im, (mmax_im/mmax_re));
  }
  
}

/*
void
write_small(double *f)
{
  FILE *fp;
  int i, j, k;
  char fname[1024];
  double mmax, mmin, maxtmp, mintmp;
  maxtmp = -1e30;
  mintmp = 1e30;

  for (k = 0; k < NZ; k++) {
    for (j = 0; j < NY; j++) {
      for (i = 0; i < NX; i++) {
        double tmp = f[k * NX * NY + j * NX + i];
        if (maxtmp < tmp) {
          maxtmp = tmp;
        }
        if (mintmp > tmp) {
          mintmp = tmp;
        }
      }
    }
  }

  MPI_Allreduce (&maxtmp, &mmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&mintmp, &mmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  // z section
  sprintf(fname, "./Bn/Bnz%d%d%d", cart_id[0], cart_id[1], cart_id[2]);
  //if (cart_id[2] == (procs[2]/2)) {
  if (cart_id[2] == 0)
  {
    fp = fopen(fname, "wb");
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
        short tmp = (short)((f[k * NX * NY + j * NX + i] - mmin) / (mmax - mmin) * 65535);
        fwrite (&tmp, sizeof (short), 1, fp);
      }
    }
    fclose(fp);
  }
  // x section
  sprintf(fname, "./Bn/Bnx%d%d%d", cart_id[1], cart_id[2], cart_id[0]);
  if (cart_id[0] == (procs[0] / 2))
  {
    fp = fopen(fname, "wb");
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
        short tmp = (short)((f[k * NX * NY + j * NX + i] - mmin) / (mmax - mmin) * 65535);
        fwrite (&tmp, sizeof (short), 1, fp);
      }
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
        short tmp = (short)((f[k * NX * NY + j * NX + i] - mmin) / (mmax - mmin) * 65535);
        fwrite (&tmp, sizeof (short), 1, fp);
      }
    }
    fclose(fp);
  }
}*/

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

void
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
void
elastic_malloc()
{
  int p, q;
  for (p = 0; p < NELASTIC; p++)
  {
  	for (q = 0; q < NELASTIC; q++)
        {
        	//Bn[p][q] = (double *)_mm_malloc(sizeof(double) * NX * NY * NZ,256);
        	Bn1[p][q] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
                memset(Bn1[p][q], 0.0, Approx * Approx * Approx * 8 * sizeof(double));
        }
  	tmpy_fftre[p] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
  	tmpy_fftim[p] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
  	Bn_fftre[p] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
  	Bn_fftim[p] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
	fftre[p] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
	fftim[p] = (double *)_mm_malloc(sizeof(double) * Approx * Approx * Approx * 8, 256);
  }
  //n0 = (double *)_mm_malloc(sizeof(double) * Approx, 256);
  //omega_cpu = (double *)_mm_malloc(sizeof(double) * DIM * DIM, 256);
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
  for (p = 0; p < NELASTIC; p++)
  {
 	for (q = 0; q < NELASTIC; q++)
	{
		//hipMalloc (&BN[p][q], NX * NY * NZ * sizeof(double));
		hipMalloc (&BN1[p][q], Approx * Approx * Approx * 8 * sizeof(double));
   //             hipMalloc (&tmp0[p][q], DIM * DIM * sizeof(double));
   //             hipMalloc (&tmp1[p][q], DIM * DIM * sizeof(double));
   //             hipMalloc (&tmp2[p][q], DIM * DIM * sizeof(double));
  //              hipMemset(tmp0[p][q], 0.0, DIM * DIM * sizeof(double));
  //              hipMemset(tmp1[p][q], 0.0, DIM * DIM * sizeof(double));
  //              hipMemset(tmp2[p][q], 0.0, DIM * DIM * sizeof(double));
	}
  //	hipMalloc (&Epsilon2d[p], DIM * DIM * sizeof(double));
  //	hipMalloc (&Sigma2d[p], DIM * DIM * sizeof(double));
//	hipMemset(Sigma2d[p], 0.0, sizeof(double) * DIM * DIM);

  	hipMalloc (&tmpy_fftRE[p], Approx * Approx * Approx * 8 * sizeof(double));
  	hipMalloc (&tmpy_fftIM[p], Approx * Approx * Approx * 8 * sizeof(double));
  	hipMalloc (&BN_fftRE[p], Approx * Approx * Approx * 8 * sizeof(double));
  	hipMalloc (&BN_fftIM[p], Approx * Approx * Approx * 8 * sizeof(double));
	hipMalloc (&fftRE[p], Approx * Approx * Approx * 8 * sizeof(double));
  	hipMalloc (&fftIM[p], Approx * Approx * Approx * 8 * sizeof(double));
  }
// hipMalloc (&n11, Approx * sizeof(double));
// hipMalloc (&n22, Approx * sizeof(double));
// hipMalloc (&n33, Approx * sizeof(double));
// hipMalloc (&n123, DIM * sizeof(double));
// hipMalloc (&omega_inver, DIM * DIM * sizeof(double));
// hipMalloc (&E, DIM * DIM * sizeof(double));
// hipMalloc (&omega, DIM * DIM * sizeof(double));
// hipMalloc (&C4D, DIM * DIM * DIM * DIM * sizeof(double));
  hipMalloc (&tmpy_RE1, NX * NY * NZ*sizeof(double));
  hipMalloc (&Gx_re, NX * Approx  * sizeof(double));
  hipMalloc (&Gx_im, NX * Approx  * sizeof(double));
  hipMalloc (&Gy_re, NX * Approx  * sizeof(double));
  hipMalloc (&Gy_im, NY * Approx  * sizeof(double));
  hipMalloc (&Gz_re, NZ * Approx  * sizeof(double));
  hipMalloc (&Gz_im, NZ * Approx  * sizeof(double));
  hipMalloc (&temp1, Approx * NY * NZ * sizeof(double));
  hipMalloc (&temp2, Approx * NY * NZ * sizeof(double));
  hipMalloc (&temp3, Approx * Approx * NZ * sizeof(double));
  hipMalloc (&temp4, Approx * Approx * NZ * sizeof(double));
  hipMalloc (&temp5, Approx * Approx * NZ * sizeof(double));
  hipMalloc (&temp6, Approx * Approx * NZ * sizeof(double));
  hipMalloc (&temp7, Approx * Approx * Approx * sizeof(double));
  hipMalloc (&temp8, Approx * Approx * Approx * sizeof(double));
  hipMalloc (&Itemp1, Approx * NY * NZ  * sizeof(double));
  hipMalloc (&Itemp2, Approx * NY * NZ  * sizeof(double));
  hipMalloc (&Itemp3, NX * NY * NZ * sizeof(double));
  hipMalloc (&Itemp4, NX * NY * NZ * sizeof(double));
  //hipMalloc (&Itemp5, NX * NY * NZ * sizeof(double));
  hipMalloc (&Elas, nx * ny * nz* sizeof (double));
  hipMalloc (&S_left, (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&R_left, (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&S_right, (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&R_right, (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&S_top, nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&R_top, nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&S_bottom, nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&R_bottom, nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof (double));
  hipMalloc (&S_front, nx * ny * (nghost + 2) * sizeof (double));
  hipMalloc (&R_front, nx * ny * (nghost + 2) * sizeof (double));
  hipMalloc (&S_back, nx * ny * (nghost + 2) * sizeof (double));
  hipMalloc (&R_back, nx * ny * (nghost + 2) * sizeof (double));
}

void
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
  MPI_Send_init(&s_front[0], 1, front_back, front, 9, MPI_COMM_WORLD, &ireq_front_back[0]);
  MPI_Recv_init(&r_back[0], 1, front_back, back, 9, MPI_COMM_WORLD, &ireq_front_back[1]);
  MPI_Send_init(&s_back[0], 1, front_back, back, 9, MPI_COMM_WORLD, &ireq_front_back[2]);
  MPI_Recv_init(&r_front[0], 1, front_back, front, 9, MPI_COMM_WORLD, &ireq_front_back[3]);
}

void elastic_finish()
{
  int p, q, l;
//  fft_finish();
  for (p = 0; p < NELASTIC; p++)
  {
    for (q = 0; q < NELASTIC; q++)
    {
      //_mm_free(Bn[p][q]);
      _mm_free(Bn1[p][q]);
      //hipFree(BN[p][q]);
      hipFree(BN1[p][q]);
 //    hipFree(tmp0[p][q]);
 //    hipFree(tmp1[p][q]);
 //    hipFree(tmp2[p][q]);
    }
    _mm_free(tmpy_fftre[p]);
    _mm_free(tmpy_fftim[p]);
    _mm_free(Bn_fftre[p]);
    _mm_free(Bn_fftim[p]);
    _mm_free(fftre[p]);
    _mm_free(fftim[p]);

//   hipFree(Epsilon2d[p]);
//   hipFree(Sigma2d[p]);

    hipFree(tmpy_fftRE[p]);
    hipFree(tmpy_fftIM[p]);
    hipFree(BN_fftRE[p]);
    hipFree(BN_fftIM[p]);
    hipFree(fftRE[p]);
    hipFree(fftIM[p]);
  }
  _mm_free(tmpy_re);
  _mm_free(tmpy_im);
// _mm_free(n0);
// _mm_free(omega_cpu);
// hipFree(C4D);
// hipFree(n11);
// hipFree(n22);
// hipFree(n33);
// hipFree(n123);
// hipFree(omega_inver);
// hipFree(omega);
// hipFree(E);

  hipFree(tmpy_RE1);
  hipFree(Gx_re);
  hipFree(Gx_im);
  hipFree(Gy_re);
  hipFree(Gy_im);
  hipFree(Gz_re);
  hipFree(Gz_im);
  hipFree(temp1);
  hipFree(temp2);
  hipFree(temp3);
  hipFree(temp4);
  hipFree(temp5);
  hipFree(temp6);
  hipFree(temp7);
  hipFree(temp8);
  hipFree(Itemp1);
  hipFree(Itemp2);
  hipFree(Itemp3);
  hipFree(Itemp4);
  //hipFree(Itemp5);
  /* free elastic variable*/
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
  hipFree(Elas);
  hipFree(S_left);
  hipFree(R_left);
  hipFree(S_right);
  hipFree(R_right);
  hipFree(S_top);
  hipFree(R_top);
  hipFree(S_bottom);
  hipFree(R_bottom);
  hipFree(S_front);
  hipFree(R_front);
  hipFree(S_back);
  hipFree(R_back);
}

void
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
  char filename[1024];
  FILE *fp;
  char *aline = NULL;
  size_t len = 0;
  double C2d[6 * 6];
  NX = nx - 2 * nghost;
  NY = ny - 2 * nghost;
  NZ = nz - 2 * nghost;
  sprintf(filename, "elastic_input.txt");
  fp = fopen(filename, "r");
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
  //hipMemcpy (C4D, C4d, DIM * DIM * DIM * DIM * sizeof (double), hipMemcpyHostToDevice);
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
    //hipMemcpy (Epsilon2d[p], epsilon2d[p], DIM * DIM * sizeof (double), hipMemcpyHostToDevice);
    // sigma0[ij] = C[ijkl] * epsilon0[kl]
   #if 1
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
    //write_n(sigma2d[p], p);
    //hipMemcpy (Sigma2d[p], sigma2d[p], DIM * DIM * sizeof (double), hipMemcpyHostToDevice);
    #endif
  }
  fclose(fp);
  if (aline)
    free(aline);
}

// omega = 1/normal(omega_inver) * (-1)^(i+j) * M(omega_inver)
void
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
void
elastic_calculate_BN()
{
  int p, q;
  int x, y, z;
  int i, j, k, l;
  double n11_local[2 * Approx], n22_local[2 * Approx], n33_local[2 * Approx];
  double normal;
  int gx, gy, gz;
  int num = Approx;

  // calculate n
  int gnnx = procs[0] * NX; // 256
  int gnny = procs[1] * NY; // 256
  int gnnz = procs[2] * NZ; // 256
  int cntx = gnnx / 2;      // 128
  int cnty = gnny / 2;      // 128
  int cntz = gnnz / 2;      // 128
  for (x = 0; x < Approx; x++)
  {
    gx = cart_id[0] * NX + x;
    n11_local[x] = 1.0 * gx;
  }
  for (x = 0; x < Approx; x++)
  {
    gx = cart_id[0] * NX + x + (NX - num);
    n11_local[x + num] = -1.0 * (gnnx - gx);
  }
  for (y = 0; y < Approx; y++)
  {
    gy = cart_id[1] * NY + y;
    n22_local[y] = 1.0 * gy;
  }
  for (y = 0; y < Approx; y++)
  {
    gy = cart_id[1] * NY + y + (NY - num);
    n22_local[y + num] = -1.0 * (gnny - gy);
  }
  for (z = 0; z < Approx; z++)
  {
    gz = cart_id[2] * NZ + z;
    n33_local[z] = 1.0 * gz;
  }
  for (z = 0; z < Approx; z++)
  {
    gz = cart_id[2] * NZ + z + (NZ - num);
    n33_local[z + num] = -1.0 * (gnnz - gz);
  }
#if 1

  /* 
    Bn(p,q) = C[ijkl] * epsilon0(p)[ij] * epsilon0(q)[kl] 
     - n[i] * sigma0(p)[ij] * omega[jk] * sigma0(q)[kl] * n[l]
  */
  double n123_local[DIM];
  double omega_inver_local[DIM * DIM];
  double omega_local[DIM * DIM];
  double tmp0_local, tmp1_local;
  double I;
  for (int stride = 0; stride < 8; stride++)
  {
    for (z = 0; z < Approx; z++)
    {
      for (y = 0; y < Approx; y++)
      {
        for (x = 0; x < Approx; x++)
        {
 
          gx = cart_id[0] * NX + x + (num * (stride > 3));
          gy = cart_id[1] * NY + y + (num * (stride == 2 || stride == 3 || stride == 6 || stride == 7));
          gz = cart_id[2] * NZ + z + (num * (stride % 2 == 1));
          if (gx == 0 && gy == 0 && gz == 0)
          {
            n123_local[0] = 0;
            n123_local[1] = 1;
            n123_local[2] = 0;
          }
          else
          {
            I = 1.0 /pow((n11_local[x+(num*(stride>3))]*n11_local[x+(num*(stride>3))]) + (n22_local[y+(num*(stride==2||stride==3||stride==6||stride==7))]*n22_local[y+(num*(stride==2||stride==3||stride==6||stride==7))]) + (n33_local[z+(num*(stride%2==1))]*n33_local[z+(num*(stride%2==1))]), 0.5);
            n123_local[0] = n11_local[x+(num*(stride>3))] * I;
            n123_local[1] = n22_local[y+(num*(stride==2||stride==3||stride==6||stride==7))] * I;
            n123_local[2] = n33_local[z+(num*(stride%2==1))] * I;
          }
 
          // omega_inver[ik] = C[ijkl] * n[j] * n[l]
          for (i = 0; i < DIM; ++i)
          {
            for (k = 0; k < DIM; ++k)
            {
              omega_inver_local[i * DIM + k] = 0.0;
              for (j = 0; j < DIM; ++j)
              {
                for (l = 0; l < DIM; ++l)
                {
                  omega_inver_local[i * DIM + k] += C4d[i * DIM * DIM * DIM + j * DIM * DIM + k * DIM + l] * n123_local[j] * n123_local[l];
                }
              }
            }
          }
 
          // omega = inverse(omega_inver)
          inverse_3x3(omega_local, omega_inver_local);
          for (p = 0; p < NELASTIC; p++)
          {
            for (q = 0; q < NELASTIC; q++)
            {
              tmp0_local = 0;
              tmp1_local = 0;
              // tmp0 = C[ijkl] * epsilon[ij] * epsilon[kl]
              // tmp1 = n[i] * sigma0[ij] * omega[jk] * sigma0[kl] *n[l]
              for (i = 0; i < DIM; i++)
              {
                for (j = 0; j < DIM; j++)
                {
                  for (k = 0; k < DIM; k++)
                  {
                    for (l = 0; l < DIM; l++)
                    {
                      tmp0_local += C4d[i * DIM * DIM * DIM + j * DIM * DIM + k * DIM + l] * epsilon2d[p][i * DIM + j] * epsilon2d[q][k * DIM + l];
                      tmp1_local += n123_local[i] * sigma2d[p][i * DIM + j] * omega_local[j * DIM + k] * sigma2d[q][k * DIM + l] * n123_local[l];
                    }
                  }
                }
              }
              // Bn = tmp0 - tmp1
              Bn1[p][q][z * Approx * Approx + y * Approx + x + stride * ela_size] = tmp0_local - tmp1_local;
            }
      	  }
        }
      }
    }
  }
  for (p = 0; p < NELASTIC; p++)
  {
    for (q = 0; q < NELASTIC; q++)
    {
      hipMemcpy (BN1[p][q], Bn1[p][q], Approx * Approx * Approx * 8 * sizeof (double), hipMemcpyHostToDevice);
    }
  }
#endif
}
#if 0
void
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
//      n11[x] = 2.0*PI/gnnx/hx*gx; // 0-127
    }
    else
    {                              // 128-255
      n11[x] = -1.0 * (gnnx - gx); // (-128)-(-1)
//      n11[x] = 2.0*PI/gnnx/hx*(gx-gnnx); // (-128)-(-1)
    }
    //printf("n11[%d] = %+1.15lf\n", x, n11[x]);
  }

  for (y = 0; y < NY; y++) // 0-64
  {
    int gy = cart_id[1] * NY + y; // 0-255
    if (gy < cnty)                // 0-127
    {
      n22[y] = 1.0 * gy; // 0-127
//      n22[y] = 2.0*PI/gnny/hy*gy;
    }
    else
    {                              // 128-255
      n22[y] = -1.0 * (gnny - gy); // (-128)-(-1)
//      n22[y] = 2.0*PI/gnny/hy*(gy-gnny); // (-128)-(-1)
    }
  }
  for (z = 0; z < NZ; z++) // 0-64
  {
    int gz = cart_id[2] * NZ + z; // 0-255
    if (gz < cntz)                // 0-127
    {
      n33[z] = 1.0 * gz; // 0-127
//      n33[z] = 2.0*PI/gnnz/hz*gz; // 0-127
    }
    else
    {                              // 128-255
      n33[z] = -1.0 * (gnnz - gz); // (-128)-(-1)
//      n33[z] = 2.0*PI/gnnz/hz*(gz-gnnz); // (-128)-(-1)
    }
  }
  /* 
    Bn(p,q) = C[ijkl] * epsilon0(p)[ij] * epsilon0(q)[kl] 
     - n[i] * sigma0(p)[ij] * omega[jk] * sigma0(q)[kl] * n[l]
  */
  for (p = 0; p < 12; p++)
  {
    for (q = 0; q < 12; q++)
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
	//write_n(omega_inver, p, q);

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
            //if (gx == 0 && gy == 0 && gz == 0)
	    //{
            //Bn[p][q][z * NY * NX + y * NX + x] = tmp0;
	    //}
	    //else
	    //{
            Bn[p][q][z * NY * NX + y * NX + x] = tmp0 - tmp1;
	    //}
	    //Bn[z * NY * NX + y * NX + x] = 75 - tmp1;
          }
        }
      }
	write_Bn_re(Bn[p][q], p, q);
//	Bn[p][q][0] = 0.0;
//      hipMemcpy (BN[p][q], Bn[p][q], NZ* NY* NX* sizeof (double), hipMemcpyHostToDevice);
    }
  }
  for (p = 1; p < NELASTIC; p++)
  {
    for (q = 1; q < NELASTIC; q++)
    {
//      hipMemcpy (BN[p][q], Bn[p][q], NZ* NY* NX* sizeof (double), hipMemcpyHostToDevice);
    }
  }
}
#endif
void elastic_init()
{
  //fft_setup();
  elastic_malloc();
  elastic_init_transfer();
  elastic_calculate_BN();
}

