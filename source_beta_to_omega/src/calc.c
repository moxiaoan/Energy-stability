#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "ScLETD.h"
#include "mkl.h"


void ac_calc_FU(int n, double *f)
{
  int i, j, k;
  // ac[n].f1
  ac_calc_F1(n);
  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        f[k * ny * nx + j * nx + i] = ac[n].f1[k * ny * nx + j * nx + i];
      }
    }
  }
  // ac[n].f2
  if (ANISOTROPIC == 1)
  {
    ac_calc_F2(n);
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        for (i = ix1; i < ix4; i++)
        {
          f[k * ny * nx + j * nx + i] -= epn2 * ac[n].LE * ac[n].f2[k * ny * nx + j * nx + i];
          //epn2=0.6955
        }
      }
    }
  }
  // ac[n].elas
  if (ELASTIC == AC_FUNCTION)
  {
//    printf("elastic\n");
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        for (i = ix1; i < ix4; i++)
        {
          f[k * ny * nx + j * nx + i] -= ac[n].LE * ElasticScale * ac[n].felas[k * ny * nx + j * nx + i];
        }
      }
    }
  }
  // random noise: ac[n].f3
  if (iter < 300)
  {
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        for (i = ix1; i < ix4; i++)
        {
          //f[k * ny * nx + j * nx + i] += ac[n].f3[k * ny * nx + j * nx + i];
        }
      }
    }
  }
}

void ch_calc_FU_ConstantMobility(int n, double *fieldci1)
{
  int m;
  int i, j, k;
  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double f_left, f_right, f_top, f_bottom, f_front, f_back;
  double tmp1, tmp2;

  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      //#pragma simd
      for (i = ix1; i < ix4; i++)
      {
        SWITCH_CH_FIELD1(ch[n].fieldCIt[k * nx * ny + j * nx + i]);
      }
    }
  }

  if (ELASTIC == CH_FUNCTION)
  {
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        for (i = ix1; i < ix4; i++)
        {
          ch[n].fieldCIt[k * ny * nx + j * nx + i] += ch[n].LCI * ch[n].felas[k * ny * nx + j * nx + i];
        }
      }
    }
  }

  if (left >= 0)
  {
    i = ix1;
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        ch[n].fieldCIt[k * nx * ny + j * nx + i] -= ch[n].LCI * epn2 * ch[n].fieldCIu_left[k * ny + j] / hx / hx;
      }
    }
  }

  if (right >= 0)
  {
    i = ix4 - 1;
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        ch[n].fieldCIt[k * nx * ny + j * nx + i] -= ch[n].LCI * epn2 * ch[n].fieldCIu_right[k * ny + j] / hx / hx;
      }
    }
  }

  if (top >= 0)
  {
    j = iy1;
    for (k = iz1; k < iz4; k++)
    {
      for (i = ix1; i < ix4; i++)
      {
        ch[n].fieldCIt[k * nx * ny + j * nx + i] -= ch[n].LCI * epn2 * ch[n].fieldCIu_top[k * nx + i] / hy / hy;
      }
    }
  }

  if (bottom >= 0)
  {
    j = iy4 - 1;
    for (k = iz1; k < iz4; k++)
    {
      for (i = ix1; i < ix4; i++)
      {
        ch[n].fieldCIt[k * nx * ny + j * nx + i] -= ch[n].LCI * epn2 * ch[n].fieldCIu_bottom[k * nx + i] / hy / hy;
      }
    }
  }

  if (front >= 0)
  {
    k = iz1;
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        ch[n].fieldCIt[k * nx * ny + j * nx + i] -= ch[n].LCI * epn2 * ch[n].fieldCIu_front[j * nx + i] / hz / hz;
      }
    }
  }

  if (back >= 0)
  {
    k = iz4 - 1;
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        ch[n].fieldCIt[k * nx * ny + j * nx + i] -= ch[n].LCI * epn2 * ch[n].fieldCIu_back[j * nx + i] / hz / hz;
      }
    }
  }

  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        if (left < 0)
        {
          if (i == ix1)
          {
            n_left = ix1 + 1;
            n_right = ix1 + 1;
          }
          else if (i < ix4 - 1)
          {
            n_left = i - 1;
            n_right = i + 1;
          }
          else
          {
            n_left = i - 1;
            n_right = -1;
          }
        }
        else if (right >= 0)
        {
          if (i == ix1)
          {
            n_left = -1;
            n_right = i + 1;
          }
          else if (i < ix4 - 1)
          {
            n_left = i - 1;
            n_right = i + 1;
          }
          else
          {
            n_left = i - 1;
            n_right = -1;
          }
        }
        else
        {
          if (i == ix1)
          {
            n_left = -1;
            n_right = i + 1;
          }
          else if (i < ix4 - 1)
          {
            n_left = i - 1;
            n_right = i + 1;
          }
          else
          {
            n_left = i - 1;
            n_right = i - 1;
          }
        }

        if (top < 0)
        {
          if (j == iy1)
          {
            n_top = j + 1;
            n_bottom = j + 1;
          }
          else if (j < iy4 - 1)
          {
            n_top = j - 1;
            n_bottom = j + 1;
          }
          else
          {
            n_top = j - 1;
            n_bottom = -1;
          }
        }
        else if (bottom >= 0)
        {
          if (j == iy1)
          {
            n_top = -1;
            n_bottom = j + 1;
          }
          else if (j < iy4 - 1)
          {
            n_top = j - 1;
            n_bottom = j + 1;
          }
          else
          {
            n_top = j - 1;
            n_bottom = -1;
          }
        }
        else
        {
          if (j == iy1)
          {
            n_top = -1;
            n_bottom = j + 1;
          }
          else if (j < iy4 - 1)
          {
            n_top = j - 1;
            n_bottom = j + 1;
          }
          else
          {
            n_top = j - 1;
            n_bottom = j - 1;
          }
        }
        if (front < 0)
        {
          if (k == iz1)
          {
            n_front = k + 1;
            n_back = k + 1;
          }
          else if (k < iz4 - 1)
          {
            n_front = k - 1;
            n_back = k + 1;
          }
          else
          {
            n_front = k - 1;
            n_back = -1;
          }
        }
        else if (back >= 0)
        {
          if (k == iz1)
          {
            n_front = -1;
            n_back = k + 1;
          }
          else if (k < iz4 - 1)
          {
            n_front = k - 1;
            n_back = k + 1;
          }
          else
          {
            n_front = k - 1;
            n_back = -1;
          }
        }
        else
        {
          if (k == iz1)
          {
            n_front = -1;
            n_back = k + 1;
          }
          else if (k < iz4 - 1)
          {
            n_front = k - 1;
            n_back = k + 1;
          }
          else
          {
            n_front = k - 1;
            n_back = k - 1;
          }
        }

        if (n_front > -1)
        {
          f_front = ch[n].fieldCIt[n_front * nx * ny + j * nx + i];
        }
        else
        {
          f_front = 0;
        }
        if (n_back > -1)
        {
          f_back = ch[n].fieldCIt[n_back * nx * ny + j * nx + i];
        }
        else
        {
          f_back = 0;
        }
        if (n_top > -1)
        {
          f_top = ch[n].fieldCIt[k * nx * ny + n_top * nx + i];
        }
        else
        {
          f_top = 0;
        }
        if (n_bottom > -1)
        {
          f_bottom = ch[n].fieldCIt[k * nx * ny + n_bottom * nx + i];
        }
        else
        {
          f_bottom = 0;
        }
        if (n_right > -1)
        {
          f_right = ch[n].fieldCIt[k * nx * ny + j * nx + n_right];
        }
        else
        {
          f_right = 0;
        }
        if (n_left > -1)
        {
          f_left = ch[n].fieldCIt[k * nx * ny + j * nx + n_left];
        }
        else
        {
          f_left = 0;
        }

        fieldci1[k * nx * ny + j * nx + i] = (f_left + f_right - 2.0 * ch[n].fieldCIt[k * nx * ny + j * nx + i]) / hx / hx;
        fieldci1[k * nx * ny + j * nx + i] += (f_top + f_bottom - 2.0 * ch[n].fieldCIt[k * nx * ny + j * nx + i]) / hy / hy;
        fieldci1[k * nx * ny + j * nx + i] += (f_front + f_back - 2.0 * ch[n].fieldCIt[k * nx * ny + j * nx + i]) / hz / hz;
      }
    }
  }

  if (left >= 0)
  {
    i = ix1;
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        fieldci1[k * nx * ny + j * nx + i] -= ch[n].fieldCImu_left[k * ny + j] / hx / hx;
      }
    }
  }

  if (right >= 0)
  {
    i = ix4 - 1;
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        fieldci1[k * nx * ny + j * nx + i] -= ch[n].fieldCImu_right[k * ny + j] / hx / hx;
      }
    }
  }

  if (top >= 0)
  {
    j = iy1;
    for (k = iz1; k < iz4; k++)
    {
      for (i = ix1; i < ix4; i++)
      {
        fieldci1[k * nx * ny + j * nx + i] -= ch[n].fieldCImu_top[k * nx + i] / hy / hy;
      }
    }
  }

  if (bottom >= 0)
  {
    j = iy4 - 1;
    for (k = iz1; k < iz4; k++)
    {
      for (i = ix1; i < ix4; i++)
      {
        fieldci1[k * nx * ny + j * nx + i] -= ch[n].fieldCImu_bottom[k * nx + i] / hy / hy;
      }
    }
  }

  if (front >= 0)
  {
    k = iz1;
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        fieldci1[k * nx * ny + j * nx + i] -= ch[n].fieldCImu_front[j * nx + i] / hz / hz;
      }
    }
  }

  if (back >= 0)
  {
    k = iz4 - 1;
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        fieldci1[k * nx * ny + j * nx + i] -= ch[n].fieldCImu_back[j * nx + i] / hz / hz;
      }
    }
  }
}

void ch_calc_FU(int n, double *f)
{
  if (VariableMobility)
  {
    update_M_C();
    ch_calc_FU_VariableMobility(n, f);
  }
  else
  {
    printf("ch_calc_FU_ConstantMobility\n");
    ch_calc_FU_ConstantMobility(n, f);
  }
}

void ac_updateU_new(int n, double *field, double *field1)
{
  double tmp, Hijk;
  int i, j, k, l;
  for (j = iy1; j < iy4; j++)
  {
    for (i = ix1; i < ix4; i++)
    {
      for (k = iz1; k < iz4; k++)
      {
        l = j * nz * nx + i * nz + k;
        tmp = kkz * DDZ[k] + kky * DDY[j] + kkx * DDX[i];
        Hijk = -ac[n].LE * (tmp * epn2 - ac[n].KE);
        if (fabs(Hijk) < 1.0e-8)
        {
          Hijk = 0.0;
        }
        tmp = 1.0 - ac[n].phiE[l] * Hijk;
        field[l] = tmp * field[l] + ac[n].phiE[l] * field1[l];
      }
    }
  }
}

void ch_updateU_new(int n, double *field, double *field1)
{
  double tmp, Hijk;
  int i, j, k, l;
  for (j = iy1; j < iy4; j++)
  {
    for (i = ix1; i < ix4; i++)
    {
      for (k = iz1; k < iz4; k++)
      {
        l = j * nz * nx + i * nz + k;
        tmp = DDZ[k] + DDY[j] + DDX[i];
        Hijk = ch[n].LCI * (tmp * tmp * epn2 - ch[n].KCI * tmp);
        if (fabs(Hijk) < 1.0e-8)
        {
          Hijk = 0.0;
        }
        tmp = 1.0 - ch[n].phiCI[l] * Hijk;
        field[l] = tmp * field[l] + ch[n].phiCI[l] * field1[l];
      }
    }
  }
}

prepare_U1_new(double *field1, double *field2)
{
  int i, j, k, l;
  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        l = k * nx * ny + j * nx + i;
        field2[l] = field2[l] - field1[l];
      }
    }
  }
}

void prepare_U2_new(double *phi, double *field1, double *field2)
{
  int i, j, k, l;
  for (j = iy1; j < iy4; j++)
  {
    for (i = ix1; i < ix4; i++)
    {
      for (k = iz1; k < iz4; k++)
      {
        l = j * nx * nz + i * nz + k;
        field2[l] = phi[l] * field2[l];
      }
    }
  }
}

void correct_U_new(double *field, double *field1)
{
  int i, j, k, l;
  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        l = k * nx * ny + j * nx + i;
        field[l] += field1[l];
      }
    }
  }
}

void xyz_yzx(double *f, double *ft)
{
  int i, j, k;
  for (i = ix1; i < ix4; i++)
  {
    for (k = iz1; k < iz4; k++)
    {
      for (j = iy1; j < iy4; j++)
      {
        ft[i * ny * nz + k * ny + j] = f[k * nx * ny + j * nx + i];
      }
    }
  }
}

void yzx_zxy(double *f, double *ft)
{
  int i, j, k;
  for (j = iy1; j < iy4; j++)
  {
    for (i = ix1; i < ix4; i++)
    {
      for (k = iz1; k < iz4; k++)
      {
        ft[j * nz * nx + i * nz + k] = f[i * ny * nz + k * ny + j];
      }
    }
  }
}

void zxy_xyz(double *f, double *ft)
{
  int i, j, k;

  for (k = iz1; k < iz4; k++)
  {
    for (j = iy1; j < iy4; j++)
    {
      for (i = ix1; i < ix4; i++)
      {
        ft[k * nx * ny + j * nx + i] = f[j * nz * nx + i * nz + k];
      }
    }
  }
}

void PUX(double *A, double *B, double *C, double *D)
{
  int m, n, k;

  m = nx;
  n = ny * nz;
  k = nx;

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, m, B, k, beta, D, m);
  xyz_yzx(D, C);
}

void PUY(double *A, double *B, double *C, double *D)
{
  int m, n, k;

  m = ny;
  n = nz * nx;
  k = ny;

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, m, B, k, beta, D, m);
  yzx_zxy(D, C);
}

void PUZ(double *A, double *B, double *C, double *D, double *E)
{
  int m, n, k;

  m = nz;
  n = nx * ny;
  k = nz;

  switch (stage)
  {
  case 0:

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, m, B, k, beta, C, m);

    break;

  case 1:

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, m, B, k, beta, D, m);

    break;

  case 2:

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, m, B, k, beta, D, m);
    xyz_yzx(D, C);

    break;
  }
}
