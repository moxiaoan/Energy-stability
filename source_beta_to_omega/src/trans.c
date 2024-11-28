#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "ScLETD.h"
#include "mkl.h"

void top_bottom_pack(double *field, double *fields_top, double *fields_bottom, double *fieldr_front, double *fieldr_back)
{
  int k_s, k_e;
  int f_j, f_k;
  int tb_j, tb_k;
  int t_ofst, b_ofst, i, j, k;

  k_s = nghost + 2;
  k_e = nz + (nghost + 2);
  f_j = nx;
  f_k = nx * ny;
  tb_j = nx;
  tb_k = nx * (nghost + 2);
  t_ofst = nx * nghost;
  b_ofst = nx * (ny - nghost - (nghost + 2));

  for (k = 0; k < nz + (nghost + 2) * 2; k++)
  {
    if (k < k_s)
    {
      if (front >= 0)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 0; i < nx; i++)
          {
            fields_top[tb_k * k + tb_j * j + i] = fieldr_front[f_k * k + f_j * j + i + t_ofst];
            fields_bottom[tb_k * k + tb_j * j + i] = fieldr_front[f_k * k + f_j * j + i + b_ofst];
          }
        }
      }
    }
    else if (k >= k_e)
    {
      if (back >= 0)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 0; i < nx; i++)
          {
            fields_top[tb_k * k + tb_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * j + i + t_ofst];
            fields_bottom[tb_k * k + tb_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * j + i + b_ofst];
          }
        }
      }
    }
    else
    {
      for (j = 0; j < nghost + 2; j++)
      {
        for (i = 0; i < nx; i++)
        {
          fields_top[tb_k * k + tb_j * j + i] = field[f_k * (k - k_s) + f_j * j + i + t_ofst];
          fields_bottom[tb_k * k + tb_j * j + i] = field[f_k * (k - k_s) + f_j * j + i + b_ofst];
        }
      }
    }
  }
}

void left_right_pack(double *field, double *fields_left, double *fields_right, double *fieldr_top, double *fieldr_bottom, double *fieldr_front, double *fieldr_back)
{
  int j_s, j_e, k_s, k_e;
  int f_j, f_k;
  int lr_j, lr_k;
  int tb_j, tb_k;
  int l_ofst, r_ofst, i, j, k;

  j_s = nghost + 2;
  j_e = ny + (nghost + 2);
  k_s = nghost + 2;
  k_e = nz + (nghost + 2);
  f_j = nx;
  f_k = nx * ny;
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);
  tb_j = nx;
  tb_k = nx * (nghost + 2);
  l_ofst = nghost;
  r_ofst = nx - nghost - (nghost + 2);

  for (k = 0; k < nz + (nghost + 2) * 2; k++)
  {
    if (top >= 0)
    {
      for (j = 0; j < j_s; j++)
      {
        for (i = 0; i < nghost + 2; i++)
        {
          fields_left[lr_k * k + lr_j * j + i] = fieldr_top[tb_k * k + tb_j * j + i + l_ofst];
          fields_right[lr_k * k + lr_j * j + i] = fieldr_top[tb_k * k + tb_j * j + i + r_ofst];
        }
      }
    }
    if (bottom >= 0)
    {
      for (j = j_e; j < ny + (nghost + 2) * 2; j++)
      {
        for (i = 0; i < nghost + 2; i++)
        {
          fields_left[lr_k * k + lr_j * j + i] = fieldr_bottom[tb_k * k + tb_j * (j - j_e) + i + l_ofst];
          fields_right[lr_k * k + lr_j * j + i] = fieldr_bottom[tb_k * k + tb_j * (j - j_e) + i + r_ofst];
        }
      }
    }
    if (k < k_s)
    {
      if (front >= 0)
      {
        for (j = j_s; j < j_e; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fields_left[lr_k * k + lr_j * j + i] = fieldr_front[f_k * k + f_j * (j - j_s) + i + l_ofst];
            fields_right[lr_k * k + lr_j * j + i] = fieldr_front[f_k * k + f_j * (j - j_s) + i + r_ofst];
          }
        }
      }
    }
    else if (k >= k_e)
    {
      if (back >= 0)
      {
        for (j = j_s; j < j_e; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fields_left[lr_k * k + lr_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * (j - j_s) + i + l_ofst];
            fields_right[lr_k * k + lr_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * (j - j_s) + i + r_ofst];
          }
        }
      }
    }
    else
    {
      for (j = j_s; j < j_e; j++)
      {
        for (i = 0; i < nghost + 2; i++)
        {
          fields_left[lr_k * k + lr_j * j + i] = field[f_k * (k - k_s) + f_j * (j - j_s) + i + l_ofst];
          fields_right[lr_k * k + lr_j * j + i] = field[f_k * (k - k_s) + f_j * (j - j_s) + i + r_ofst];
        }
      }
    }
  }
}

void unpack(double *field, double *fieldr_left, double *fieldr_right, double *fieldr_top, double *fieldr_bottom, double *fieldr_front, double *fieldr_back)
{
  int f_j, f_k;
  int lr_k, lr_j;
  int tb_k, tb_j;
  int i, j, k;

  f_j = nx;
  f_k = nx * ny;
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);
  tb_j = nx;
  tb_k = nx * (nghost + 2);

  if (left >= 0)
  {
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nghost; i++)
        {
          field[f_k * k + f_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i + 2];
        }
      }
    }
  }

  if (right >= 0)
  {
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nghost; i++)
        {
          field[f_k * k + f_j * j + i + (nx - nghost)] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i];
        }
      }
    }
  }

  if (top >= 0)
  {
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < nghost; j++)
      {
        for (i = 0; i < nx; i++)
        {
          field[f_k * k + f_j * j + i] = fieldr_top[tb_k * (k + (nghost + 2)) + tb_j * (j + 2) + i];
        }
      }
    }
  }

  if (bottom >= 0)
  {
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < nghost; j++)
      {
        for (i = 0; i < nx; i++)
        {
          field[f_k * k + f_j * (j + ny - nghost) + i] = fieldr_bottom[tb_k * (k + (nghost + 2)) + tb_j * j + i];
        }
      }
    }
  }

  if (front >= 0)
  {
    for (k = 0; k < nghost; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          field[f_k * k + f_j * j + i] = fieldr_front[f_k * (k + 2) + f_j * j + i];
        }
      }
    }
  }

  if (back >= 0)
  {
    for (k = 0; k < nghost; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          field[f_k * (k + nz - nghost) + f_j * j + i] = fieldr_back[f_k * k + f_j * j + i];
        }
      }
    }
  }
}

void enlarge(double *fielde_left, double *fielde_right, double *fielde_top, double *fielde_bottom, double *fielde_front, double *fielde_back,
             double *fieldr_left, double *fieldr_right, double *fieldr_top, double *fieldr_bottom, double *fieldr_front, double *fieldr_back)
{
  int j_s, j_e, k_s, k_e;
  int elr_k, elr_j, etb_k, etb_j, efb_k, efb_j;
  int lr_k, lr_j, tb_k, tb_j, fb_k, fb_j;
  int b_ofst;
  int i, j, k;

  j_s = nghost + 2;
  j_e = j_s + ny - (nghost + 2);
  if (top < 0)
  {
    j_s -= nghost;
  }
  if (bottom < 0)
  {
    j_e += nghost;
  }
  k_s = nghost + 2;
  k_e = k_s + nz - (nghost + 2);
  if (front < 0)
  {
    k_s -= nghost;
  }
  if (back < 0)
  {
    k_e += nghost;
  }

  elr_j = nghost + 2;
  elr_k = (nghost + 2) * (ny + 4);
  etb_j = nx + 4;
  etb_k = (nx + 4) * (nghost + 2);
  efb_j = (nx + 4) * (nghost + 2);
  efb_k = nx + 4;
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);
  tb_j = nx;
  tb_k = nx * (nghost + 2);
  fb_j = nx;
  fb_k = nx * ny;
  b_ofst = nz + (nghost + 2);

  if (left >= 0)
  {
    for (k = 0; k < nz + 4; k++)
    {
      if (front >= 0 && k < k_s)
      {
        for (j = 2; j < ny + 2; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + nghost) + i];
          }
        }

        if (top >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * k + lr_j * j + i];
            }
          }
        }

        if (bottom >= 0)
        {
          for (j = ny; j < ny + 4; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }
      }

      if (k >= k_s && k < k_e)
      {
        for (j = 2; j < ny + 2; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * (j + nghost) + i];
          }
        }

        if (top >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * j + i];
            }
          }
        }

        if (bottom >= 0)
        {
          for (j = ny; j < ny + 4; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }
      }

      if (back >= 0 && k >= k_e)
      {
        for (j = 2; j < ny + 2; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + nghost) + i];
          }
        }

        if (top >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * j + i];
            }
          }
        }

        if (bottom >= 0)
        {
          for (j = ny; j < ny + 4; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }
      }
    }
  }
  if (right >= 0)
  {
    for (k = 0; k < nz + 4; k++)
    {
      if (front >= 0 && k < k_s)
      {
        for (j = 2; j < ny + 2; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + nghost) + i];
          }
        }

        if (top >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * k + lr_j * j + i];
            }
          }
        }

        if (bottom >= 0)
        {
          for (j = ny; j < ny + 4; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }
      }

      if (k >= k_s && k < k_e)
      {
        for (j = 2; j < ny + 2; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * (j + nghost) + i];
          }
        }

        if (top >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * j + i];
            }
          }
        }

        if (bottom >= 0)
        {
          for (j = ny; j < ny + 4; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }
      }

      if (back >= 0 && k >= k_e)
      {
        for (j = 2; j < ny + 2; j++)
        {
          for (i = 0; i < nghost + 2; i++)
          {
            fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + nghost) + i];
          }
        }

        if (top >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * j + i];
            }
          }
        }

        if (bottom >= 0)
        {
          for (j = ny; j < ny + 4; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }
      }
    }
  }

  if (top >= 0)
  {
    for (k = 0; k < nz + 4; k++)
    {
      if (front >= 0 && k < k_s)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_top[etb_k * k + etb_j * j + i] = fieldr_top[tb_k * k + tb_j * j + (i - 2)];
          }
        }

        if (left >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_top[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * k + lr_j * j + i];
            }
          }
        }

        if (right >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_top[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * k + lr_j * j + (i - nx)];
            }
          }
        }
      }

      if (k >= k_s && k < k_e)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_top[etb_k * k + etb_j * j + i] = fieldr_top[tb_k * (k + nghost) + tb_j * j + (i - 2)];
          }
        }

        if (left >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_top[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * j + i];
            }
          }
        }

        if (right >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_top[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * j + (i - nx)];
            }
          }
        }
      }

      if (back >= 0 && k >= k_e)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_top[etb_k * k + etb_j * j + i] = fieldr_top[tb_k * (k + (nghost + 2)) + tb_j * j + (i - 2)];
          }
        }

        if (left >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_top[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * j + i];
            }
          }
        }

        if (right >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_top[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * j + (i - nx)];
            }
          }
        }
      }
    }
  }

  if (bottom >= 0)
  {
    for (k = 0; k < nz + 4; k++)
    {
      if (front >= 0 && k < nghost + 2)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_bottom[tb_k * k + tb_j * j + (i - 2)];
          }
        }

        if (left >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + (ny + (nghost + 2))) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + (ny + (nghost + 2))) + (i - nx)];
            }
          }
        }
      }

      if (k >= k_s && k < k_e)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_bottom[tb_k * (k + nghost) + tb_j * j + (i - 2)];
          }
        }

        if (left >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * (j + (ny + (nghost + 2))) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * (j + (ny + (nghost + 2))) + (i - nx)];
            }
          }
        }
      }

      if (back >= 0 && k >= k_e)
      {
        for (j = 0; j < nghost + 2; j++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_bottom[tb_k * (k + (nghost + 2)) + tb_j * j + (i - 2)];
          }
        }

        if (left >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + (ny + (nghost + 2))) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (j = 0; j < nghost + 2; j++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + (ny + (nghost + 2))) + (i - nx)];
            }
          }
        }
      }
    }
  }

  if (front >= 0)
  {
    for (j = 0; j < ny + 4; j++)
    {
      if (top >= 0 && j < j_s)
      {
        for (k = 0; k < nghost + 2; k++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_front[efb_k * k + efb_j * j + i] = fieldr_top[tb_k * k + tb_j * j + i - 2];
          }
        }

        if (left >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_front[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * k + lr_j * j + i];
            }
          }
        }

        if (right >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_front[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * k + lr_j * j + i - nx];
            }
          }
        }
      }

      if (j >= j_s && j < j_e)
      {
        for (k = 0; k < nghost + 2; k++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_front[efb_k * k + efb_j * j + i] = fieldr_front[fb_k * k + fb_j * (j - nghost) + i - 2];
          }
        }

        if (left >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_front[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + nghost) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_front[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + nghost) + i - nx];
            }
          }
        }
      }

      if (bottom >= 0 && j >= j_e)
      {
        for (k = 0; k < nghost + 2; k++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_front[efb_k * k + efb_j * j + i] = fieldr_bottom[tb_k * k + tb_j * (j - j_e) + i - 2];
          }
        }

        if (left >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_front[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_front[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + (nghost + 2)) + i - nx];
            }
          }
        }
      }
    }
  }

  if (back >= 0)
  {
    for (j = 0; j < ny + 4; j++)
    {
      if (top >= 0 && j < j_s)
      {
        for (k = 0; k < nghost + 2; k++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_back[efb_k * k + efb_j * j + i] = fieldr_top[tb_k * (k + b_ofst) + tb_j * j + i - 2];
          }
        }

        if (left >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_back[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * (k + b_ofst) + lr_j * j + i];
            }
          }
        }

        if (right >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_back[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * (k + b_ofst) + lr_j * j + i - nx];
            }
          }
        }
      }

      if (j >= j_s && j < j_e)
      {
        for (k = 0; k < nghost + 2; k++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_back[efb_k * k + efb_j * j + i] = fieldr_back[fb_k * k + fb_j * (j - nghost) + i - 2];
          }
        }

        if (left >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_back[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * (k + b_ofst) + lr_j * (j + nghost) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_back[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * (k + b_ofst) + lr_j * (j + nghost) + i - nx];
            }
          }
        }
      }

      if (bottom >= 0 && j >= j_e)
      {
        for (k = 0; k < nghost + 2; k++)
        {
          for (i = 2; i < nx + 2; i++)
          {
            fielde_back[efb_k * k + efb_j * j + i] = fieldr_bottom[tb_k * (k + b_ofst) + tb_j * (j - j_e) + i - 2];
          }
        }

        if (left >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = 0; i < nghost + 2; i++)
            {
              fielde_back[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * (k + b_ofst) + lr_j * (j + (nghost + 2)) + i];
            }
          }
        }

        if (right >= 0)
        {
          for (k = 0; k < nghost + 2; k++)
          {
            for (i = nx; i < nx + 4; i++)
            {
              fielde_back[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * (k + b_ofst) + lr_j * (j + (nghost + 2)) + i - nx];
            }
          }
        }
      }
    }
  }
}

void ac_mu(double *fieldmu_left, double *fieldmu_right, double *fieldmu_top, double *fieldmu_bottom, double *fieldmu_front, double *fieldmu_back,
           double *fielde_left, double *fielde_right, double *fielde_top, double *fielde_bottom, double *fielde_front, double *fielde_back,
           double *fieldu_left, double *fieldu_right, double *fieldu_top, double *fieldu_bottom, double *fieldu_front, double *fieldu_back)
{
  int i, j, k;
  int mu_i, mu_j, mu_k;
  int e_i, e_j, e_k;
  int l_e, l_mu;

  mu_k = ny;
  e_j = nghost + 2;
  e_k = (nghost + 2) * (ny + 4);

  if (left >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      for (j = 2; j < ny + 2; j++)
      {
        l_e = k * e_k + j * e_j + 1;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_left[l_mu] = (fielde_left[l_e - 1] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + 1]) / hx / hx;
      }

      if (top < 0)
      {
        l_e = k * e_k + 2 * e_j + 1;
        l_mu = (k - 2) * mu_k;
        fieldmu_left[l_mu] += ((-2.0) * fielde_left[l_e] + (2.0) * fielde_left[l_e + e_j]) / hy / hy;
        for (j = 3; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_j]) / hy / hy;
        }
      }
      else if (bottom < 0)
      {
        for (j = 2; j < ny + 1; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_j]) / hy / hy;
        }
        l_e = k * e_k + (ny + 1) * e_j + 1;
        l_mu = (k - 2) * mu_k + (ny + 1) - 2;
        fieldmu_left[l_mu] += ((2.0) * fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e]) / hy / hy;
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_j]) / hy / hy;
        }
      }

      if (k == 2 && front < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += ((-2.0) * fielde_left[l_e] + (2.0) * fielde_left[l_e + e_k]) / hz / hz;
          fieldmu_left[l_mu] = epn2 * fieldmu_left[l_mu] - fielde_left[l_e] * (fielde_left[l_e] * fielde_left[l_e] - 1.0);
          fieldu_left[l_mu] = fielde_left[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += ((2.0) * fielde_left[l_e - e_k] + (-2.0) * fielde_left[l_e]) / hz / hz;
          fieldmu_left[l_mu] = epn2 * fieldmu_left[l_mu] - fielde_left[l_e] * (fielde_left[l_e] * fielde_left[l_e] - 1.0);
          fieldu_left[l_mu] = fielde_left[l_e];
        }
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_k] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_k]) / hz / hz;
          fieldmu_left[l_mu] = epn2 * fieldmu_left[l_mu] - fielde_left[l_e] * (fielde_left[l_e] * fielde_left[l_e] - 1.0);
          fieldu_left[l_mu] = fielde_left[l_e];
        }
      }
    }
  }

  if (right >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      for (j = 2; j < ny + 2; j++)
      {
        l_e = k * e_k + j * e_j + 2;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_right[l_mu] = (fielde_right[l_e - 1] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + 1]) / hx / hx;
      }

      if (top < 0)
      {
        l_e = k * e_k + 2 * e_j + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_right[l_mu] += ((-2.0) * fielde_right[l_e] + (2.0) * fielde_right[l_e + e_j]) / hy / hy;
        for (j = 3; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_j]) / hy / hy;
        }
      }
      else if (bottom < 0)
      {
        for (j = 2; j < ny + 1; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_j]) / hy / hy;
        }
        l_e = k * e_k + (ny + 1) * e_j + 2;
        l_mu = (k - 2) * mu_k + (ny + 1) - 2;
        fieldmu_right[l_mu] += ((2.0) * fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e]) / hy / hy;
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_j]) / hy / hy;
        }
      }

      if (k == 2 && front < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += ((-2.0) * fielde_right[l_e] + (2.0) * fielde_right[l_e + e_k]) / hz / hz;
          fieldmu_right[l_mu] = epn2 * fieldmu_right[l_mu] - fielde_right[l_e] * (fielde_right[l_e] * fielde_right[l_e] - 1.0);
          fieldu_right[l_mu] = fielde_right[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += ((2.0) * fielde_right[l_e - e_k] + (-2.0) * fielde_right[l_e]) / hz / hz;
          fieldmu_right[l_mu] = epn2 * fieldmu_right[l_mu] - fielde_right[l_e] * (fielde_right[l_e] * fielde_right[l_e] - 1.0);
          fieldu_right[l_mu] = fielde_right[l_e];
        }
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_k] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_k]) / hz / hz;
          fieldmu_right[l_mu] = epn2 * fieldmu_right[l_mu] - fielde_right[l_e] * (fielde_right[l_e] * fielde_right[l_e] - 1.0);
          fieldu_right[l_mu] = fielde_right[l_e];
        }
      }
    }
  }

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  if (top >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      if (left < 0)
      {
        l_e = k * e_k + e_j + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_top[l_mu] = ((-2.0) * fielde_top[l_e] + (2.0) * fielde_top[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] = (fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] = (fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_top[l_mu] = ((2.0) * fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] = (fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + 1]) / hx / hx;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_top[l_mu] += (fielde_top[l_e - e_j] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + e_j]) / hy / hy;
      }

      if (k == 2 && front < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] += ((-2.0) * fielde_top[l_e] + (2.0) * fielde_top[l_e + e_k]) / hz / hz;
          fieldmu_top[l_mu] = epn2 * fieldmu_top[l_mu] - fielde_top[l_e] * (fielde_top[l_e] * fielde_top[l_e] - 1.0);
          fieldu_top[l_mu] = fielde_top[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] += ((2.0) * fielde_top[l_e - e_k] + (-2.0) * fielde_top[l_e]) / hz / hz;
          fieldmu_top[l_mu] = epn2 * fieldmu_top[l_mu] - fielde_top[l_e] * (fielde_top[l_e] * fielde_top[l_e] - 1.0);
          fieldu_top[l_mu] = fielde_top[l_e];
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] += (fielde_top[l_e - e_k] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + e_k]) / hz / hz;
          fieldmu_top[l_mu] = epn2 * fieldmu_top[l_mu] - fielde_top[l_e] * (fielde_top[l_e] * fielde_top[l_e] - 1.0);
          fieldu_top[l_mu] = fielde_top[l_e];
        }
      }
    }
  }

  if (bottom >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      if (left < 0)
      {
        l_e = k * e_k + e_j * nghost + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_bottom[l_mu] = ((-2.0) * fielde_bottom[l_e] + (2.0) * fielde_bottom[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] = (fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] = (fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j * nghost + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_bottom[l_mu] = ((2.0) * fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] = (fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + 1]) / hx / hx;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j * nghost + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_bottom[l_mu] += (fielde_bottom[l_e - e_j] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + e_j]) / hy / hy;
      }

      if (k == 2 && front < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] += ((-2.0) * fielde_bottom[l_e] + (2.0) * fielde_bottom[l_e + e_k]) / hz / hz;
          fieldmu_bottom[l_mu] = epn2 * fieldmu_bottom[l_mu] - fielde_bottom[l_e] * (fielde_bottom[l_e] * fielde_bottom[l_e] - 1.0);
          fieldu_bottom[l_mu] = fielde_bottom[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] += ((2.0) * fielde_bottom[l_e - e_k] + (-2.0) * fielde_bottom[l_e]) / hz / hz;
          fieldmu_bottom[l_mu] = epn2 * fieldmu_bottom[l_mu] - fielde_bottom[l_e] * (fielde_bottom[l_e] * fielde_bottom[l_e] - 1.0);
          fieldu_bottom[l_mu] = fielde_bottom[l_e];
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] += (fielde_bottom[l_e - e_k] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + e_k]) / hz / hz;
          fieldmu_bottom[l_mu] = epn2 * fieldmu_bottom[l_mu] - fielde_bottom[l_e] * (fielde_bottom[l_e] * fielde_bottom[l_e] - 1.0);
          fieldu_bottom[l_mu] = fielde_bottom[l_e];
        }
      }
    }
  }

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  if (front >= 0)
  {
    for (k = 2; k < ny + 2; k++)
    {

      if (left < 0)
      {
        l_e = k * e_k + e_j + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_front[l_mu] = ((-2.0) * fielde_front[l_e] + (2.0) * fielde_front[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] = (fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] = (fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_front[l_mu] = ((2.0) * fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] = (fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + 1]) / hx / hx;
        }
      }

      if (k == 2 && top < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] += ((-2.0) * fielde_front[l_e] + (2.0) * fielde_front[l_e + e_k]) / hy / hy;
        }
      }
      else if (k == ny + 1 && bottom < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] += ((2.0) * fielde_front[l_e - e_k] + (-2.0) * fielde_front[l_e]) / hy / hy;
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] += (fielde_front[l_e - e_k] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + e_k]) / hy / hy;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_front[l_mu] += (fielde_front[l_e - e_j] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + e_j]) / hz / hz;
        fieldmu_front[l_mu] = epn2 * fieldmu_front[l_mu] - fielde_front[l_e] * (fielde_front[l_e] * fielde_front[l_e] - 1.0);
        fieldu_front[l_mu] = fielde_front[l_e];
      }
    }
  }

  if (back >= 0)
  {
    for (k = 2; k < ny + 2; k++)
    {

      if (left < 0)
      {
        l_e = k * e_k + e_j * nghost + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_back[l_mu] = ((-2.0) * fielde_back[l_e] + (2.0) * fielde_back[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] = (fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] = (fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j * nghost + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_back[l_mu] = ((2.0) * fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] = (fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + 1]) / hx / hx;
        }
      }

      if (k == 2 && top < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] += ((-2.0) * fielde_back[l_e] + (2.0) * fielde_back[l_e + e_k]) / hy / hy;
        }
      }
      else if (k == ny + 1 && bottom < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] += ((2.0) * fielde_back[l_e - e_k] + (-2.0) * fielde_back[l_e]) / hy / hy;
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] += (fielde_back[l_e - e_k] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + e_k]) / hy / hy;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j * nghost + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_back[l_mu] += (fielde_back[l_e - e_j] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + e_j]) / hz / hz;
        fieldmu_back[l_mu] = epn2 * fieldmu_back[l_mu] - fielde_back[l_e] * (fielde_back[l_e] * fielde_back[l_e] - 1.0);
        fieldu_back[l_mu] = fielde_back[l_e];
      }
    }
  }
}

void ch_mu(int n, double LCI, double *fieldmu_left, double *fieldmu_right, double *fieldmu_top, double *fieldmu_bottom, double *fieldmu_front, double *fieldmu_back,
           double *fielde_left, double *fielde_right, double *fielde_top, double *fielde_bottom, double *fielde_front, double *fielde_back,
           double *fieldu_left, double *fieldu_right, double *fieldu_top, double *fieldu_bottom, double *fieldu_front, double *fieldu_back)
{
  int m;
  int i, j, k;
  int mu_i, mu_j, mu_k;
  int e_i, e_j, e_k;
  int l_e, l_mu;

  mu_k = ny;
  e_j = nghost + 2;
  e_k = (nghost + 2) * (ny + 4);

  if (left >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      for (j = 2; j < ny + 2; j++)
      {
        l_e = k * e_k + j * e_j + 1;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_left[l_mu] = (fielde_left[l_e - 1] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + 1]) / hx / hx;
      }

      if (top < 0)
      {
        l_e = k * e_k + 2 * e_j + 1;
        l_mu = (k - 2) * mu_k;
        fieldmu_left[l_mu] += ((-2.0) * fielde_left[l_e] + (2.0) * fielde_left[l_e + e_j]) / hy / hy;
        for (j = 3; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_j]) / hy / hy;
        }
      }
      else if (bottom < 0)
      {
        for (j = 2; j < ny + 1; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_j]) / hy / hy;
        }
        l_e = k * e_k + (ny + 1) * e_j + 1;
        l_mu = (k - 2) * mu_k + (ny + 1) - 2;
        fieldmu_left[l_mu] += ((2.0) * fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e]) / hy / hy;
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_j] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_j]) / hy / hy;
        }
      }

      if (k == 2 && front < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += ((-2.0) * fielde_left[l_e] + (2.0) * fielde_left[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(left);
          fieldu_left[l_mu] = fielde_left[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += ((2.0) * fielde_left[l_e - e_k] + (-2.0) * fielde_left[l_e]) / hz / hz;
          SWITCH_CH_FIELDMU(left);
          fieldu_left[l_mu] = fielde_left[l_e];
        }
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_left[l_mu] += (fielde_left[l_e - e_k] + (-2.0) * fielde_left[l_e] + fielde_left[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(left);
          fieldu_left[l_mu] = fielde_left[l_e];
        }
      }
    }
  }

  if (right >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      for (j = 2; j < ny + 2; j++)
      {
        l_e = k * e_k + j * e_j + 2;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_right[l_mu] = (fielde_right[l_e - 1] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + 1]) / hx / hx;
      }

      if (top < 0)
      {
        l_e = k * e_k + 2 * e_j + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_right[l_mu] += ((-2.0) * fielde_right[l_e] + (2.0) * fielde_right[l_e + e_j]) / hy / hy;
        for (j = 3; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_j]) / hy / hy;
        }
      }
      else if (bottom < 0)
      {
        for (j = 2; j < ny + 1; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_j]) / hy / hy;
        }
        l_e = k * e_k + (ny + 1) * e_j + 2;
        l_mu = (k - 2) * mu_k + (ny + 1) - 2;
        fieldmu_right[l_mu] += ((2.0) * fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e]) / hy / hy;
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_j] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_j]) / hy / hy;
        }
      }

      // front and back
      if (k == 2 && front < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += ((-2.0) * fielde_right[l_e] + (2.0) * fielde_right[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(right);
          fieldu_right[l_mu] = fielde_right[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += ((2.0) * fielde_right[l_e - e_k] + (-2.0) * fielde_right[l_e]) / hz / hz;
          SWITCH_CH_FIELDMU(right);
          fieldu_right[l_mu] = fielde_right[l_e];
        }
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_right[l_mu] += (fielde_right[l_e - e_k] + (-2.0) * fielde_right[l_e] + fielde_right[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(right);
          //fieldmu_right[l_mu] = epn2 * fieldmu_right[l_mu] - fielde_right[l_e] * (fielde_right[l_e] * fielde_right[l_e] - 1.0);
          //fieldmu_right[l_mu] = LCI * epn2 * fieldmu_right[l_mu] - LCI * (A1 * (fielde_right[l_e] - C1) + 0.5 * A5 * ac[0].fieldEe_right[l_e] * ac[0].fieldEe_right[l_e]);
          fieldu_right[l_mu] = fielde_right[l_e];
        }
      }
    }
  }

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  // top face
  if (top >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      // left and right
      if (left < 0)
      {
        l_e = k * e_k + e_j + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_top[l_mu] = ((-2.0) * fielde_top[l_e] + (2.0) * fielde_top[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] = (fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] = (fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_top[l_mu] = ((2.0) * fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] = (fielde_top[l_e - 1] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + 1]) / hx / hx;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_top[l_mu] += (fielde_top[l_e - e_j] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + e_j]) / hy / hy;
      }

      if (k == 2 && front < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] += ((-2.0) * fielde_top[l_e] + (2.0) * fielde_top[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(top);
          fieldu_top[l_mu] = fielde_top[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] += ((2.0) * fielde_top[l_e - e_k] + (-2.0) * fielde_top[l_e]) / hz / hz;
          SWITCH_CH_FIELDMU(top);
          fieldu_top[l_mu] = fielde_top[l_e];
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_top[l_mu] += (fielde_top[l_e - e_k] + (-2.0) * fielde_top[l_e] + fielde_top[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(top);
          fieldu_top[l_mu] = fielde_top[l_e];
        }
      }
    }
  }

  if (bottom >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      if (left < 0)
      {
        l_e = k * e_k + e_j * nghost + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_bottom[l_mu] = ((-2.0) * fielde_bottom[l_e] + (2.0) * fielde_bottom[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] = (fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] = (fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j * nghost + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_bottom[l_mu] = ((2.0) * fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] = (fielde_bottom[l_e - 1] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + 1]) / hx / hx;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j * nghost + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_bottom[l_mu] += (fielde_bottom[l_e - e_j] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + e_j]) / hy / hy;
      }
      if (k == 2 && front < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] += ((-2.0) * fielde_bottom[l_e] + (2.0) * fielde_bottom[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(bottom);
          fieldu_bottom[l_mu] = fielde_bottom[l_e];
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] += ((2.0) * fielde_bottom[l_e - e_k] + (-2.0) * fielde_bottom[l_e]) / hz / hz;
          SWITCH_CH_FIELDMU(bottom);
          fieldu_bottom[l_mu] = fielde_bottom[l_e];
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_bottom[l_mu] += (fielde_bottom[l_e - e_k] + (-2.0) * fielde_bottom[l_e] + fielde_bottom[l_e + e_k]) / hz / hz;
          SWITCH_CH_FIELDMU(bottom);
          fieldu_bottom[l_mu] = fielde_bottom[l_e];
        }
      }
    }
  }

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  if (front >= 0)
  {
    for (k = 2; k < ny + 2; k++)
    {
      if (left < 0)
      {
        l_e = k * e_k + e_j + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_front[l_mu] = ((-2.0) * fielde_front[l_e] + (2.0) * fielde_front[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] = (fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] = (fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_front[l_mu] = ((2.0) * fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] = (fielde_front[l_e - 1] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + 1]) / hx / hx;
        }
      }
      if (k == 2 && top < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] += ((-2.0) * fielde_front[l_e] + (2.0) * fielde_front[l_e + e_k]) / hy / hy;
        }
      }
      else if (k == ny + 1 && bottom < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] += ((2.0) * fielde_front[l_e - e_k] + (-2.0) * fielde_front[l_e]) / hy / hy;
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_front[l_mu] += (fielde_front[l_e - e_k] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + e_k]) / hy / hy;
        }
      }

      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_front[l_mu] += (fielde_front[l_e - e_j] + (-2.0) * fielde_front[l_e] + fielde_front[l_e + e_j]) / hz / hz;
        SWITCH_CH_FIELDMU(front);
        fieldu_front[l_mu] = fielde_front[l_e];
      }
    }
  }

  if (back >= 0)
  {
    for (k = 2; k < ny + 2; k++)
    {
      if (left < 0)
      {
        l_e = k * e_k + e_j * nghost + 2;
        l_mu = (k - 2) * mu_k;
        fieldmu_back[l_mu] = ((-2.0) * fielde_back[l_e] + (2.0) * fielde_back[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] = (fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] = (fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j * nghost + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        fieldmu_back[l_mu] = ((2.0) * fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] = (fielde_back[l_e - 1] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + 1]) / hx / hx;
        }
      }
      if (k == 2 && top < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] += ((-2.0) * fielde_back[l_e] + (2.0) * fielde_back[l_e + e_k]) / hy / hy;
        }
      }
      else if (k == ny + 1 && bottom < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] += ((2.0) * fielde_back[l_e - e_k] + (-2.0) * fielde_back[l_e]) / hy / hy;
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          fieldmu_back[l_mu] += (fielde_back[l_e - e_k] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + e_k]) / hy / hy;
        }
      }
      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j * nghost + j;
        l_mu = (k - 2) * mu_k + j - 2;
        fieldmu_back[l_mu] += (fielde_back[l_e - e_j] + (-2.0) * fielde_back[l_e] + fielde_back[l_e + e_j]) / hz / hz;
        SWITCH_CH_FIELDMU(back);
        fieldu_back[l_mu] = fielde_back[l_e];
      }
    }
  }
}

void transfer(void)
{
  int i, j, k;
  int n;
  for (n = 0; n < nac; n++)
  {
    MPI_Startall(4, ac[n].ireq_front_back_fieldE);
    MPI_Waitall(4, ac[n].ireq_front_back_fieldE, status);
  }
  for (n = 0; n < nch; n++)
  {
    MPI_Startall(4, ch[n].ireq_front_back_fieldCI);
    MPI_Waitall(4, ch[n].ireq_front_back_fieldCI, status);
  }
  for (n = 0; n < nac; n++)
  {
    top_bottom_pack(ac[n].fieldE, ac[n].fieldEs_top, ac[n].fieldEs_bottom, ac[n].fieldEr_front, ac[n].fieldEr_back);
    MPI_Startall(4, ac[n].ireq_top_bottom_fieldE);
    MPI_Waitall(4, ac[n].ireq_top_bottom_fieldE, status);
  }
  for (n = 0; n < nch; n++)
  {
    top_bottom_pack(ch[n].fieldCI, ch[n].fieldCIs_top, ch[n].fieldCIs_bottom, ch[n].fieldCIr_front, ch[n].fieldCIr_back);
    MPI_Startall(4, ch[n].ireq_top_bottom_fieldCI);
    MPI_Waitall(4, ch[n].ireq_top_bottom_fieldCI, status);
  }

  for (n = 0; n < nac; n++)
  {
    left_right_pack(ac[n].fieldE, ac[n].fieldEs_left, ac[n].fieldEs_right, ac[n].fieldEr_top, ac[n].fieldEr_bottom, ac[n].fieldEr_front, ac[n].fieldEr_back);
    MPI_Startall(4, ac[n].ireq_left_right_fieldE);
    MPI_Waitall(4, ac[n].ireq_left_right_fieldE, status);
  }
  for (n = 0; n < nch; n++)
  {
    left_right_pack(ch[n].fieldCI, ch[n].fieldCIs_left, ch[n].fieldCIs_right, ch[n].fieldCIr_top, ch[n].fieldCIr_bottom, ch[n].fieldCIr_front, ch[n].fieldCIr_back);
    MPI_Startall(4, ch[n].ireq_left_right_fieldCI);
    MPI_Waitall(4, ch[n].ireq_left_right_fieldCI, status);
  }
  for (n = 0; n < nac; n++)
  {
    unpack(ac[n].fieldE, ac[n].fieldEr_left, ac[n].fieldEr_right, ac[n].fieldEr_top, ac[n].fieldEr_bottom, ac[n].fieldEr_front, ac[n].fieldEr_back);
  }
  for (n = 0; n < nch; n++)
  {
    unpack(ch[n].fieldCI, ch[n].fieldCIr_left, ch[n].fieldCIr_right, ch[n].fieldCIr_top, ch[n].fieldCIr_bottom, ch[n].fieldCIr_front, ch[n].fieldCIr_back);
  }
  for (n = 0; n < nac; n++)
  {
    enlarge(ac[n].fieldEe_left, ac[n].fieldEe_right, ac[n].fieldEe_top, ac[n].fieldEe_bottom, ac[n].fieldEe_front, ac[n].fieldEe_back,
            ac[n].fieldEr_left, ac[n].fieldEr_right, ac[n].fieldEr_top, ac[n].fieldEr_bottom, ac[n].fieldEr_front, ac[n].fieldEr_back);
  }
  for (n = 0; n < nch; n++)
  {
    enlarge(ch[n].fieldCIe_left, ch[n].fieldCIe_right, ch[n].fieldCIe_top, ch[n].fieldCIe_bottom, ch[n].fieldCIe_front, ch[n].fieldCIe_back,
            ch[n].fieldCIr_left, ch[n].fieldCIr_right, ch[n].fieldCIr_top, ch[n].fieldCIr_bottom, ch[n].fieldCIr_front, ch[n].fieldCIr_back);
  }
}

void calc_mu()
{
  int n;
  for (n = 0; n < nac; n++)
  {
    ac_mu(ac[n].fieldEmu_left, ac[n].fieldEmu_right, ac[n].fieldEmu_top, ac[n].fieldEmu_bottom, ac[n].fieldEmu_front, ac[n].fieldEmu_back,
          ac[n].fieldEe_left, ac[n].fieldEe_right, ac[n].fieldEe_top, ac[n].fieldEe_bottom, ac[n].fieldEe_front, ac[n].fieldEe_back,
          ac[n].fieldEu_left, ac[n].fieldEu_right, ac[n].fieldEu_top, ac[n].fieldEu_bottom, ac[n].fieldEu_front, ac[n].fieldEu_back);
  }
  for (n = 0; n < nch; n++)
  {
    ch_mu(n, ch[n].LCI, ch[n].fieldCImu_left, ch[n].fieldCImu_right, ch[n].fieldCImu_top, ch[n].fieldCImu_bottom, ch[n].fieldCImu_front, ch[n].fieldCImu_back,
          ch[n].fieldCIe_left, ch[n].fieldCIe_right, ch[n].fieldCIe_top, ch[n].fieldCIe_bottom, ch[n].fieldCIe_front, ch[n].fieldCIe_back,
          ch[n].fieldCIu_left, ch[n].fieldCIu_right, ch[n].fieldCIu_top, ch[n].fieldCIu_bottom, ch[n].fieldCIu_front, ch[n].fieldCIu_back);
  }
}
