#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "ScLETD.h"

void ch_mu(int n, double L, double *Umu_left, double *Umu_right, double *Umu_top, double *Umu_bottom, double *Umu_front, double *Umu_back,
           double *Ue_left, double *Ue_right, double *Ue_top, double *Ue_bottom, double *Ue_front, double *Ue_back)
{
  int i, j, k;
  int mu_i, mu_j, mu_k;
  int e_i, e_j, e_k;
  int l_e, l_mu;

  mu_k = ny;
  e_j = nghost + 2;
  e_k = (nghost + 2) * (ny + 4);

  // left face
  if (left >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      // left and right
      for (j = 2; j < ny + 2; j++)
      {
        l_e = k * e_k + j * e_j + 1;
        l_mu = (k - 2) * mu_k + j - 2;
        Umu_left[l_mu] = (Ue_left[l_e - 1] + (-2.0) * Ue_left[l_e] + Ue_left[l_e + 1]) / hx / hx;
      }

      // top and bottom
      if (top < 0)
      {
        l_e = k * e_k + 2 * e_j + 1;
        l_mu = (k - 2) * mu_k;
        Umu_left[l_mu] += ((-2.0) * Ue_left[l_e] + (2.0) * Ue_left[l_e + e_j]) / hy / hy;
        for (j = 3; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_left[l_mu] += (Ue_left[l_e - e_j] + (-2.0) * Ue_left[l_e] + Ue_left[l_e + e_j]) / hy / hy;
        }
      }
      else if (bottom < 0)
      {
        for (j = 2; j < ny + 1; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_left[l_mu] += (Ue_left[l_e - e_j] + (-2.0) * Ue_left[l_e] + Ue_left[l_e + e_j]) / hy / hy;
        }
        l_e = k * e_k + (ny + 1) * e_j + 1;
        l_mu = (k - 2) * mu_k + (ny + 1) - 2;
        Umu_left[l_mu] += ((2.0) * Ue_left[l_e - e_j] + (-2.0) * Ue_left[l_e]) / hy / hy;
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_left[l_mu] += (Ue_left[l_e - e_j] + (-2.0) * Ue_left[l_e] + Ue_left[l_e + e_j]) / hy / hy;
        }
      }

      // front and back
      if (k == 2 && front < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_left[l_mu] += ((-2.0) * Ue_left[l_e] + (2.0) * Ue_left[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(left);
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_left[l_mu] += ((2.0) * Ue_left[l_e - e_k] + (-2.0) * Ue_left[l_e]) / hz / hz;
          SWITCH_CH_UMU(left);
        }
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 1;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_left[l_mu] += (Ue_left[l_e - e_k] + (-2.0) * Ue_left[l_e] + Ue_left[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(left);
        }
      }
    }
  }

  // right face
  if (right >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      // left and right
      for (j = 2; j < ny + 2; j++)
      {
        l_e = k * e_k + j * e_j + 2;
        l_mu = (k - 2) * mu_k + j - 2;
        Umu_right[l_mu] = (Ue_right[l_e - 1] + (-2.0) * Ue_right[l_e] + Ue_right[l_e + 1]) / hx / hx;
      }

      // top and bottom
      if (top < 0)
      {
        l_e = k * e_k + 2 * e_j + 2;
        l_mu = (k - 2) * mu_k;
        Umu_right[l_mu] += ((-2.0) * Ue_right[l_e] + (2.0) * Ue_right[l_e + e_j]) / hy / hy;
        for (j = 3; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_right[l_mu] += (Ue_right[l_e - e_j] + (-2.0) * Ue_right[l_e] + Ue_right[l_e + e_j]) / hy / hy;
        }
      }
      else if (bottom < 0)
      {
        for (j = 2; j < ny + 1; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_right[l_mu] += (Ue_right[l_e - e_j] + (-2.0) * Ue_right[l_e] + Ue_right[l_e + e_j]) / hy / hy;
        }
        l_e = k * e_k + (ny + 1) * e_j + 2;
        l_mu = (k - 2) * mu_k + (ny + 1) - 2;
        Umu_right[l_mu] += ((2.0) * Ue_right[l_e - e_j] + (-2.0) * Ue_right[l_e]) / hy / hy;
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_right[l_mu] += (Ue_right[l_e - e_j] + (-2.0) * Ue_right[l_e] + Ue_right[l_e + e_j]) / hy / hy;
        }
      }

      // front and back
      if (k == 2 && front < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_right[l_mu] += ((-2.0) * Ue_right[l_e] + (2.0) * Ue_right[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(right);
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_right[l_mu] += ((2.0) * Ue_right[l_e - e_k] + (-2.0) * Ue_right[l_e]) / hz / hz;
          SWITCH_CH_UMU(right);
        }
      }
      else
      {
        for (j = 2; j < ny + 2; j++)
        {
          l_e = k * e_k + j * e_j + 2;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_right[l_mu] += (Ue_right[l_e - e_k] + (-2.0) * Ue_right[l_e] + Ue_right[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(right);
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
        Umu_top[l_mu] = ((-2.0) * Ue_top[l_e] + (2.0) * Ue_top[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_top[l_mu] = (Ue_top[l_e - 1] + (-2.0) * Ue_top[l_e] + Ue_top[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_top[l_mu] = (Ue_top[l_e - 1] + (-2.0) * Ue_top[l_e] + Ue_top[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        Umu_top[l_mu] = ((2.0) * Ue_top[l_e - 1] + (-2.0) * Ue_top[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_top[l_mu] = (Ue_top[l_e - 1] + (-2.0) * Ue_top[l_e] + Ue_top[l_e + 1]) / hx / hx;
        }
      }

      // top and bottom
      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j + j;
        l_mu = (k - 2) * mu_k + j - 2;
        Umu_top[l_mu] += (Ue_top[l_e - e_j] + (-2.0) * Ue_top[l_e] + Ue_top[l_e + e_j]) / hy / hy;
      }

      // front and back
      if (k == 2 && front < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_top[l_mu] += ((-2.0) * Ue_top[l_e] + (2.0) * Ue_top[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(top);
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_top[l_mu] += ((2.0) * Ue_top[l_e - e_k] + (-2.0) * Ue_top[l_e]) / hz / hz;
          SWITCH_CH_UMU(top);
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_top[l_mu] += (Ue_top[l_e - e_k] + (-2.0) * Ue_top[l_e] + Ue_top[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(top);
        }
      }
    }
  }

  // bottom face
  if (bottom >= 0)
  {
    for (k = 2; k < nz + 2; k++)
    {

      // left and right
      if (left < 0)
      {
        l_e = k * e_k + e_j * nghost + 2;
        l_mu = (k - 2) * mu_k;
        Umu_bottom[l_mu] = ((-2.0) * Ue_bottom[l_e] + (2.0) * Ue_bottom[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_bottom[l_mu] = (Ue_bottom[l_e - 1] + (-2.0) * Ue_bottom[l_e] + Ue_bottom[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_bottom[l_mu] = (Ue_bottom[l_e - 1] + (-2.0) * Ue_bottom[l_e] + Ue_bottom[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j * nghost + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        Umu_bottom[l_mu] = ((2.0) * Ue_bottom[l_e - 1] + (-2.0) * Ue_bottom[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_bottom[l_mu] = (Ue_bottom[l_e - 1] + (-2.0) * Ue_bottom[l_e] + Ue_bottom[l_e + 1]) / hx / hx;
        }
      }

      // top and bottom
      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j * nghost + j;
        l_mu = (k - 2) * mu_k + j - 2;
        Umu_bottom[l_mu] += (Ue_bottom[l_e - e_j] + (-2.0) * Ue_bottom[l_e] + Ue_bottom[l_e + e_j]) / hy / hy;
      }

      // front and back
      if (k == 2 && front < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_bottom[l_mu] += ((-2.0) * Ue_bottom[l_e] + (2.0) * Ue_bottom[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(bottom);
        }
      }
      else if (k == nz + 1 && back < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_bottom[l_mu] += ((2.0) * Ue_bottom[l_e - e_k] + (-2.0) * Ue_bottom[l_e]) / hz / hz;
          SWITCH_CH_UMU(bottom);
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_bottom[l_mu] += (Ue_bottom[l_e - e_k] + (-2.0) * Ue_bottom[l_e] + Ue_bottom[l_e + e_k]) / hz / hz;
          SWITCH_CH_UMU(bottom);
        }
      }
    }
  }

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  // front face
  if (front >= 0)
  {
    for (k = 2; k < ny + 2; k++)
    {

      // left and right
      if (left < 0)
      {
        l_e = k * e_k + e_j + 2;
        l_mu = (k - 2) * mu_k;
        Umu_front[l_mu] = ((-2.0) * Ue_front[l_e] + (2.0) * Ue_front[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_front[l_mu] = (Ue_front[l_e - 1] + (-2.0) * Ue_front[l_e] + Ue_front[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_front[l_mu] = (Ue_front[l_e - 1] + (-2.0) * Ue_front[l_e] + Ue_front[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        Umu_front[l_mu] = ((2.0) * Ue_front[l_e - 1] + (-2.0) * Ue_front[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_front[l_mu] = (Ue_front[l_e - 1] + (-2.0) * Ue_front[l_e] + Ue_front[l_e + 1]) / hx / hx;
        }
      }

      // top and bottom
      if (k == 2 && top < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_front[l_mu] += ((-2.0) * Ue_front[l_e] + (2.0) * Ue_front[l_e + e_k]) / hy / hy;
        }
      }
      else if (k == ny + 1 && bottom < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_front[l_mu] += ((2.0) * Ue_front[l_e - e_k] + (-2.0) * Ue_front[l_e]) / hy / hy;
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_front[l_mu] += (Ue_front[l_e - e_k] + (-2.0) * Ue_front[l_e] + Ue_front[l_e + e_k]) / hy / hy;
        }
      }

      // front and back
      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j + j;
        l_mu = (k - 2) * mu_k + j - 2;
        Umu_front[l_mu] += (Ue_front[l_e - e_j] + (-2.0) * Ue_front[l_e] + Ue_front[l_e + e_j]) / hz / hz;
        SWITCH_CH_UMU(front);
      }
    }
  }

  // back face
  if (back >= 0)
  {
    for (k = 2; k < ny + 2; k++)
    {

      // left and right
      if (left < 0)
      {
        l_e = k * e_k + e_j * nghost + 2;
        l_mu = (k - 2) * mu_k;
        Umu_back[l_mu] = ((-2.0) * Ue_back[l_e] + (2.0) * Ue_back[l_e + 1]) / hx / hx;
        for (j = 3; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_back[l_mu] = (Ue_back[l_e - 1] + (-2.0) * Ue_back[l_e] + Ue_back[l_e + 1]) / hx / hx;
        }
      }
      else if (right < 0)
      {
        for (j = 2; j < nx + 1; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_back[l_mu] = (Ue_back[l_e - 1] + (-2.0) * Ue_back[l_e] + Ue_back[l_e + 1]) / hx / hx;
        }
        l_e = k * e_k + e_j * nghost + (nx + 1);
        l_mu = (k - 2) * mu_k + (nx + 1) - 2;
        Umu_back[l_mu] = ((2.0) * Ue_back[l_e - 1] + (-2.0) * Ue_back[l_e]) / hx / hx;
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_back[l_mu] = (Ue_back[l_e - 1] + (-2.0) * Ue_back[l_e] + Ue_back[l_e + 1]) / hx / hx;
        }
      }

      // top and bottom
      if (k == 2 && top < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_back[l_mu] += ((-2.0) * Ue_back[l_e] + (2.0) * Ue_back[l_e + e_k]) / hy / hy;
        }
      }
      else if (k == ny + 1 && bottom < 0)
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_back[l_mu] += ((2.0) * Ue_back[l_e - e_k] + (-2.0) * Ue_back[l_e]) / hy / hy;
        }
      }
      else
      {
        for (j = 2; j < nx + 2; j++)
        {
          l_e = k * e_k + e_j * nghost + j;
          l_mu = (k - 2) * mu_k + j - 2;
          Umu_back[l_mu] += (Ue_back[l_e - e_k] + (-2.0) * Ue_back[l_e] + Ue_back[l_e + e_k]) / hy / hy;
        }
      }

      // front and back
      for (j = 2; j < nx + 2; j++)
      {
        l_e = k * e_k + e_j * nghost + j;
        l_mu = (k - 2) * mu_k + j - 2;
        Umu_back[l_mu] += (Ue_back[l_e - e_j] + (-2.0) * Ue_back[l_e] + Ue_back[l_e + e_j]) / hz / hz;
        SWITCH_CH_UMU(back);
      }
    }
  }
}

