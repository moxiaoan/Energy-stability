#include "ScLETD.h"
#include "anisotropic_hip.h"
#include <assert.h>
#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"
#include "rocfft.h"
#include "hipfft.h"
#include <complex.h>
#include <rocfft-export.h>
#include <rocfft-version.h>
#define SCLETD_DEBUG
#define Vm 10e-10
#define al 1.07e11

__global__ void
anisotropic_calc_dev_rc(double *f, double *fa, double *lambda, int n)
{
  int i, j, k;
  //  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double a_left, a_right, a_top, a_bottom, a_front, a_back, a_middle;
  double a_left_top, a_left_bottom, a_left_front, a_left_back;
  double a_right_top, a_right_bottom, a_right_front, a_right_back;
  double a_top_front, a_top_back, a_bottom_front, a_bottom_back;
  double fxx, fyy, fzz, fxy, fxz, fyz;
  double a_top_temp, a_bottom_temp, a_front_temp, a_back_temp;
  double a_top_front_temp, a_top_back_temp;
  double a_bottom_front_temp, a_bottom_back_temp;

  int f_k = (ny + 2 * 2) * (nx + 2 * 2);
  int f_j = (nx + 2 * 2);
  int offset = 2 * f_k + 2 * f_j + 2;

  // Id of thread in the warp.
  int WARPSIZE_X = hipBlockDim_x;
  int WARPSIZE_Y = hipBlockDim_y;
  int WARPSIZE_Z = hipBlockDim_z;
  int localId_x = hipThreadIdx_x % WARPSIZE_X;
  int localId_y = hipThreadIdx_y % WARPSIZE_Y;
  int localId_z = hipThreadIdx_z % WARPSIZE_Z;

  // The first index of output element computed by this warp.
  int startOfWarp_x = hipBlockDim_x * hipBlockIdx_x + WARPSIZE_X * (hipThreadIdx_x / WARPSIZE_X);
  int startOfWarp_y = hipBlockDim_y * hipBlockIdx_y + WARPSIZE_Y * (hipThreadIdx_y / WARPSIZE_Y);
  int startOfWarp_z = hipBlockDim_z * hipBlockIdx_z + WARPSIZE_Z * (hipThreadIdx_z / WARPSIZE_Z);

  // The Id of the thread in the scope of the grid.
  i = localId_x + startOfWarp_x;
  j = localId_y + startOfWarp_y;
  k = localId_z + startOfWarp_z;
  // int globalId = k * mesh.nx * mesh.ny + j * mesh.nx + i;

  //////////////--left--middle--right///////////////////////////////////////////////////
  double rc[19];
  int globalId = k * f_k + j * f_j + i + offset;
  rc[0] = fa[k * f_k + j * f_j + i + offset];

  // if (localId_x ==(WARPSIZE_X-1) && globalId - WARPSIZE_X >=0 ) //left
  if (localId_x == (WARPSIZE_X - 1) && globalId - offset >= 0)
  {
    rc[1] = fa[globalId - WARPSIZE_X];                        // left
    rc[7] = fa[globalId - (nx + 4) - WARPSIZE_X];             // left_top
    rc[9] = fa[globalId + (nx + 4) - WARPSIZE_X];             // left_bottom
    rc[11] = fa[globalId - (ny + 4) * (nx + 4) - WARPSIZE_X]; // left_front
    rc[13] = fa[globalId + (ny + 4) * (nx + 4) - WARPSIZE_X]; // left_back
  }
  if (localId_x == 0 && WARPSIZE_X + globalId < (nx + 4) * (ny + 4) * (nz + 4))
  {
    rc[2] = fa[WARPSIZE_X + globalId];                        // right
    rc[8] = fa[globalId - (nx + 4) + WARPSIZE_X];             // right_top
    rc[10] = fa[globalId + (nx + 4) + WARPSIZE_X];            // right_bottom
    rc[12] = fa[globalId - (ny + 4) * (nx + 4) + WARPSIZE_X]; // right_front
    rc[14] = fa[globalId + (ny + 4) * (nx + 4) + WARPSIZE_X]; // right_back
  }
  if (localId_y == (WARPSIZE_Y - 1) && globalId - offset >= 0)
  {
    rc[3] = fa[globalId - WARPSIZE_Y * (nx + 4)];                        // top
    rc[15] = fa[globalId - (ny + 4) * (nx + 4) - WARPSIZE_Y * (nx + 4)]; // top_front
    rc[16] = fa[globalId + (ny + 4) * (nx + 4) - WARPSIZE_Y * (nx + 4)]; // top_back
  }
  if (localId_y == 0 && globalId + WARPSIZE_Y * (nx + 4) < (nx + 4) * (ny + 4) * (nz + 4))
  {
    rc[4] = fa[globalId + WARPSIZE_Y * (nx + 4)];                        // bottom
    rc[17] = fa[globalId - (ny + 4) * (nx + 4) + WARPSIZE_Y * (nx + 4)]; // bottom_back
    rc[18] = fa[globalId + (ny + 4) * (nx + 4) + WARPSIZE_Y * (nx + 4)]; // bottom_back
  }
  if (localId_z == (WARPSIZE_Z - 1) && globalId - offset >= 0)
  {
    rc[5] = fa[globalId - WARPSIZE_Z * (ny + 4) * (nx + 4)]; // front
  }
  if (localId_z == 0 && globalId + WARPSIZE_Z * (ny + 4) * (nx + 4) < (nx + 4) * (ny + 4) * (nz + 4))
  {
    rc[6] = fa[WARPSIZE_Z * (ny + 4) * (nx + 4) + globalId]; // back
  }
  double toShare = rc[0];
  a_middle = __shfl(toShare, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + localId_y * WARPSIZE_X + localId_x), 64);
  if (localId_x == (WARPSIZE_X - 1)) // left
  {
    toShare = rc[1];
  }
  a_left = __shfl(toShare, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + localId_y * WARPSIZE_X + localId_x - 1 + WARPSIZE_X) % WARPSIZE_X, WARPSIZE_X);
  //////
  if (localId_x == (WARPSIZE_X - 1))
  {
    toShare = rc[0];
  }
  if (localId_x == 0) // right
  {
    toShare = rc[2];
  }
  a_right = __shfl(toShare, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + localId_y * WARPSIZE_X + localId_x + 1 + WARPSIZE_X) % WARPSIZE_X, WARPSIZE_X);
  //////
  if (localId_x == 0) // right
  {
    toShare = rc[0];
  }
  if (localId_y == (WARPSIZE_Y - 1)) // top
  {
    toShare = rc[3];
  }
  a_top = __shfl(toShare, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y - 1) * WARPSIZE_X + localId_x + WARPSIZE_Y * WARPSIZE_X) % (WARPSIZE_Y * WARPSIZE_X), WARPSIZE_Y * WARPSIZE_X);
  //////
  if (localId_y == (WARPSIZE_Y - 1)) // top
  {
    toShare = rc[0];
  }
  if (localId_y == 0) // bottom
  {
    toShare = rc[4];
  }
  a_bottom = __shfl(toShare, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y + 1) * WARPSIZE_X + localId_x + WARPSIZE_Y * WARPSIZE_X) % (WARPSIZE_Y * WARPSIZE_X), WARPSIZE_Y * WARPSIZE_X);
  //////
  if (localId_y == 0) // bottom
  {
    toShare = rc[0];
  }
  if (localId_z == (WARPSIZE_Z - 1)) // front
  {
    toShare = rc[5];
  }
  a_front = __shfl(toShare, ((localId_z - 1) * WARPSIZE_Y * WARPSIZE_X + localId_y * WARPSIZE_X + localId_x), 64);
  //////
  if (localId_z == (WARPSIZE_Z - 1)) // front
  {
    toShare = rc[0];
  }
  if (localId_z == 0) // back
  {
    toShare = rc[6];
  }
  a_back = __shfl(toShare, ((localId_z + 1) * WARPSIZE_Y * WARPSIZE_X + localId_y * WARPSIZE_X + localId_x), 64);

  if (left < 0)
  {
    if (i == ix1)
    {
      a_left = a_right;
    }
  }
  if (right < 0)
  {
    if (i == ix4 - 1)
    {
      a_right = a_left;
    }
  }
  if (top < 0)
  {
    if (j == iy1)
    {
      a_top = a_bottom;
    }
  }
  if (bottom < 0)
  {
    if (j == iy4 - 1)
    {
      a_bottom = a_top;
    }
  }
  if (front < 0)
  {
    if (k == iz1)
    {
      a_front = a_back;
    }
  }
  if (back < 0)
  {
    if (k == iz4 - 1)
    {
      a_back = a_front;
    }
  }
  //////
  /*
     if (localId_z == 0) //back
     {
         toShare = rc[0]; //The following code is not used for toShare.
     }
 */
  a_top_temp = a_top;
  a_bottom_temp = a_bottom;
  a_front_temp = a_front;
  a_back_temp = a_back;
  if (localId_x == (WARPSIZE_X - 1))
  {
    a_top_temp = rc[7];    // left_top
    a_bottom_temp = rc[9]; // left_bottom
    a_front_temp = rc[11]; // left_front
    a_back_temp = rc[13];  // left_back
  }
  a_left_top = __shfl(a_top_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x - 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);
  a_left_bottom = __shfl(a_bottom_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x - 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);
  a_left_front = __shfl(a_front_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x - 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);
  a_left_back = __shfl(a_back_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x - 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);

  //////
  if (localId_x == (WARPSIZE_X - 1))
  {
    a_top_temp = a_top;       // left_top
    a_bottom_temp = a_bottom; // left_bottom
    a_front_temp = a_front;   // left_front
    a_back_temp = a_back;     // left_back
  }
  if (localId_x == 0)
  {
    a_top_temp = rc[8];     // right_top
    a_bottom_temp = rc[10]; // right_bottom
    a_front_temp = rc[12];  // right_front
    a_back_temp = rc[14];   // right_back
  }
  a_right_top = __shfl(a_top_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x + 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);
  a_right_bottom = __shfl(a_bottom_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x + 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);
  a_right_front = __shfl(a_front_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x + 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);
  a_right_back = __shfl(a_back_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y)*WARPSIZE_X + localId_x + 1 + WARPSIZE_X) % (WARPSIZE_X), WARPSIZE_X);

  //////

  a_top_front_temp = a_front;
  a_top_back_temp = a_back;
  if (localId_y == (WARPSIZE_Y - 1))
  {
    a_top_front_temp = rc[15]; // top_front
    a_top_back_temp = rc[16];  // top_back
  }
  a_top_front = __shfl(a_top_front_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y - 1) * WARPSIZE_X + localId_x + WARPSIZE_Y * WARPSIZE_X) % (WARPSIZE_Y * WARPSIZE_X), WARPSIZE_Y * WARPSIZE_X);
  a_top_back = __shfl(a_top_back_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y - 1) * WARPSIZE_X + localId_x + WARPSIZE_Y * WARPSIZE_X) % (WARPSIZE_Y * WARPSIZE_X), WARPSIZE_Y * WARPSIZE_X);
  //////
  a_bottom_front_temp = a_front;
  a_bottom_back_temp = a_back;
  if (localId_y == 0)
  {
    a_bottom_front_temp = rc[17]; // bottom_front
    a_bottom_back_temp = rc[18];  // bottom_back
  }
  a_bottom_front = __shfl(a_bottom_front_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y + 1) * WARPSIZE_X + localId_x + WARPSIZE_Y * WARPSIZE_X) % (WARPSIZE_Y * WARPSIZE_X), WARPSIZE_Y * WARPSIZE_X);
  a_bottom_back = __shfl(a_bottom_back_temp, ((localId_z)*WARPSIZE_Y * WARPSIZE_X + (localId_y + 1) * WARPSIZE_X + localId_x + WARPSIZE_Y * WARPSIZE_X) % (WARPSIZE_Y * WARPSIZE_X), WARPSIZE_Y * WARPSIZE_X);

  double tmp = 0.0;
  fxx = (a_right - 2 * a_middle + a_left) / hx / hx;
  fyy = (a_bottom - 2 * a_middle + a_top) / hy / hy;
  fzz = (a_back - 2 * a_middle + a_front) / hz / hz;
  fxy = ((a_right_bottom - a_left_bottom) - (a_right_top - a_left_top)) * 0.25 / hx / hy;
  fxz = ((a_right_back - a_left_back) - (a_right_front - a_left_front)) * 0.25 / hx / hz;
  fyz = ((a_bottom_back - a_top_back) - (a_bottom_front - a_top_front)) * 0.25 / hy / hz;
  tmp += lambda[0 * 3 + 0] * fxx;
  tmp += lambda[0 * 3 + 1] * fxy;
  tmp += lambda[0 * 3 + 2] * fxz;
  tmp += lambda[1 * 3 + 0] * fxy;
  tmp += lambda[1 * 3 + 1] * fyy;
  tmp += lambda[1 * 3 + 2] * fyz;
  tmp += lambda[2 * 3 + 0] * fxz;
  tmp += lambda[2 * 3 + 1] * fyz;
  tmp += lambda[2 * 3 + 2] * fzz;
  f[k * nx * ny + j * nx + i] = tmp;
}

__global__ void
anisotropic_calc_dev(double *f, double *fa, double *lambda, int *orie, int n)
{
  int i, j, k;
  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double a_left, a_right, a_top, a_bottom, a_front, a_back, a_middle;
  double a_left_top, a_left_bottom, a_left_front, a_left_back;
  double a_right_top, a_right_bottom, a_right_front, a_right_back;
  double a_top_front, a_top_back, a_bottom_front, a_bottom_back;
  double fxx, fyy, fzz, fxy, fxz, fyz;

  int f_k = (ny + 2 * 2) * (nx + 2 * 2);
  int f_j = (nx + 2 * 2);
  int offset = 2 * f_k + 2 * f_j + 2;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x; //  iz4, iy4, ix4 (k,j,i)
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y; //  iz1, iy1, ix1
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  n_left = i - 1;
  if (left < 0)
  {
    if (i == ix1)
    {
      n_left = i + 1;
    }
  }

  n_right = i + 1;
  if (right < 0)
  {
    if (i == ix4 - 1)
    {
      n_right = i - 1;
    }
  }

  n_top = j - 1;
  if (top < 0)
  {
    if (j == iy1)
    {
      n_top = j + 1;
    }
  }

  n_bottom = j + 1;
  if (bottom < 0)
  {
    if (j == iy4 - 1)
    {
      n_bottom = j - 1;
    }
  }

  n_front = k - 1;
  if (front < 0)
  {
    if (k == iz1)
    {
      n_front = k + 1;
    }
  }

  n_back = k + 1;
  if (back < 0)
  {
    if (k == iz4 - 1)
    {
      n_back = k - 1;
    }
  }

  a_middle = fa[k * f_k + j * f_j + i + offset];

  a_front = fa[n_front * f_k + j * f_j + i + offset];
  a_back = fa[n_back * f_k + j * f_j + i + offset];
  a_top = fa[k * f_k + n_top * f_j + i + offset];
  a_bottom = fa[k * f_k + n_bottom * f_j + i + offset];
  a_right = fa[k * f_k + j * f_j + n_right + offset];
  a_left = fa[k * f_k + j * f_j + n_left + offset];

  a_left_top = fa[k * f_k + n_top * f_j + n_left + offset];
  a_left_bottom = fa[k * f_k + n_bottom * f_j + n_left + offset];
  a_left_front = fa[n_front * f_k + j * f_j + n_left + offset];
  a_left_back = fa[n_back * f_k + j * f_j + n_left + offset];

  a_right_top = fa[k * f_k + n_top * f_j + n_right + offset];
  a_right_bottom = fa[k * f_k + n_bottom * f_j + n_right + offset];
  a_right_front = fa[n_front * f_k + j * f_j + n_right + offset];
  a_right_back = fa[n_back * f_k + j * f_j + n_right + offset];

  a_top_front = fa[n_front * f_k + n_top * f_j + i + offset];
  a_top_back = fa[n_back * f_k + n_top * f_j + i + offset];

  a_bottom_front = fa[n_front * f_k + n_bottom * f_j + i + offset];
  a_bottom_back = fa[n_back * f_k + n_bottom * f_j + i + offset];

  double tmp = 0.0;
  fxx = (a_right - 2 * a_middle + a_left) / hx / hx;
  fyy = (a_bottom - 2 * a_middle + a_top) / hy / hy;
  fzz = (a_back - 2 * a_middle + a_front) / hz / hz;
  fxy = ((a_right_bottom - a_left_bottom) - (a_right_top - a_left_top)) / 2.0 / hx / 2.0 / hy;
  fxz = ((a_right_back - a_left_back) - (a_right_front - a_left_front)) / 2.0 / hx / 2.0 / hz;
  fyz = ((a_bottom_back - a_top_back) - (a_bottom_front - a_top_front)) / 2.0 / hy / 2.0 / hz;
  int o = orie[k * nx * ny + j * nx + i];
  if (o <= 10)
  {
    tmp += lambda[0 * 3 + 0] * fxx;
    tmp += lambda[0 * 3 + 1] * fxy;
    tmp += lambda[0 * 3 + 2] * fxz;
    tmp += lambda[1 * 3 + 0] * fxy;
    tmp += lambda[1 * 3 + 1] * fyy;
    tmp += lambda[1 * 3 + 2] * fyz;
    tmp += lambda[2 * 3 + 0] * fxz;
    tmp += lambda[2 * 3 + 1] * fyz;
    tmp += lambda[2 * 3 + 2] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 10 && o <= 20)
  {
    tmp += lambda[0 * 3 + 0 + 1 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 1 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 1 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 1 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 1 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 1 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 1 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 1 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 1 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 20 && o <= 30)
  {
    tmp += lambda[0 * 3 + 0 + 3 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 3 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 3 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 3 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 3 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 3 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 3 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 3 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 3 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 30 && o <= 40)
  {
    tmp += lambda[0 * 3 + 0 + 4 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 4 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 4 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 4 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 4 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 4 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 4 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 4 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 4 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 40 && o <= 50)
  {
    tmp += lambda[0 * 3 + 0 + 5 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 5 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 5 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 5 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 5 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 5 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 5 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 5 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 5 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 50 && o <= 60)
  {
    tmp += lambda[0 * 3 + 0 + 6 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 6 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 6 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 6 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 6 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 6 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 6 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 6 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 6 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 60 && o <= 70)
  {
    tmp += lambda[0 * 3 + 0 + 7 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 7 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 7 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 7 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 7 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 7 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 7 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 7 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 7 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 70 && o <= 80)
  {
    tmp += lambda[0 * 3 + 0 + 8 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 8 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 8 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 8 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 8 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 8 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 8 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 8 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 8 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 80 && o <= 90)
  {
    tmp += lambda[0 * 3 + 0 + 9 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 9 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 9 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 9 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 9 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 9 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 9 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 9 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 9 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 90 && o <= 100)
  {
    tmp += lambda[0 * 3 + 0 + 10 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 10 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 10 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 10 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 10 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 10 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 10 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 10 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 10 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
  else if (o > 100 && o <= 125)
  {
    tmp += lambda[0 * 3 + 0 + 11 * 3 * 3] * fxx;
    tmp += lambda[0 * 3 + 1 + 11 * 3 * 3] * fxy;
    tmp += lambda[0 * 3 + 2 + 11 * 3 * 3] * fxz;
    tmp += lambda[1 * 3 + 0 + 11 * 3 * 3] * fxy;
    tmp += lambda[1 * 3 + 1 + 11 * 3 * 3] * fyy;
    tmp += lambda[1 * 3 + 2 + 11 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 0 + 11 * 3 * 3] * fxz;
    tmp += lambda[2 * 3 + 1 + 11 * 3 * 3] * fyz;
    tmp += lambda[2 * 3 + 2 + 11 * 3 * 3] * fzz;
    f[k * nx * ny + j * nx + i] = tmp;
  }
}
__global__ void
anisotropic_calc_dev_2(double *f, double *fa, double *lambda, int *orie, int n)
{
  int i, j, k;
  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double a_left, a_right, a_top, a_bottom, a_front, a_back, a_middle;
  double a_left_top, a_left_bottom, a_left_front, a_left_back;
  double a_right_top, a_right_bottom, a_right_front, a_right_back;
  double a_top_front, a_top_back, a_bottom_front, a_bottom_back;
  double fxx, fyy, fzz, fxy, fxz, fyz;

  int f_k = (ny + 2 * 2) * (nx + 2 * 2);
  int f_j = (nx + 2 * 2);
  int offset = 2 * f_k + 2 * f_j + 2;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x; //  iz4, iy4, ix4 (k,j,i)
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y; //  iz1, iy1, ix1
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  n_left = i - 1;
  if (left < 0)
  {
    if (i == ix1)
    {
      n_left = i + 1;
    }
  }

  n_right = i + 1;
  if (right < 0)
  {
    if (i == ix4 - 1)
    {
      n_right = i - 1;
    }
  }

  n_top = j - 1;
  if (top < 0)
  {
    if (j == iy1)
    {
      n_top = j + 1;
    }
  }

  n_bottom = j + 1;
  if (bottom < 0)
  {
    if (j == iy4 - 1)
    {
      n_bottom = j - 1;
    }
  }

  n_front = k - 1;
  if (front < 0)
  {
    if (k == iz1)
    {
      n_front = k + 1;
    }
  }

  n_back = k + 1;
  if (back < 0)
  {
    if (k == iz4 - 1)
    {
      n_back = k - 1;
    }
  }

  a_middle = fa[k * f_k + j * f_j + i + offset];

  a_front = fa[n_front * f_k + j * f_j + i + offset];
  a_back = fa[n_back * f_k + j * f_j + i + offset];
  a_top = fa[k * f_k + n_top * f_j + i + offset];
  a_bottom = fa[k * f_k + n_bottom * f_j + i + offset];
  a_right = fa[k * f_k + j * f_j + n_right + offset];
  a_left = fa[k * f_k + j * f_j + n_left + offset];

  a_left_top = fa[k * f_k + n_top * f_j + n_left + offset];
  a_left_bottom = fa[k * f_k + n_bottom * f_j + n_left + offset];
  a_left_front = fa[n_front * f_k + j * f_j + n_left + offset];
  a_left_back = fa[n_back * f_k + j * f_j + n_left + offset];

  a_right_top = fa[k * f_k + n_top * f_j + n_right + offset];
  a_right_bottom = fa[k * f_k + n_bottom * f_j + n_right + offset];
  a_right_front = fa[n_front * f_k + j * f_j + n_right + offset];
  a_right_back = fa[n_back * f_k + j * f_j + n_right + offset];

  a_top_front = fa[n_front * f_k + n_top * f_j + i + offset];
  a_top_back = fa[n_back * f_k + n_top * f_j + i + offset];

  a_bottom_front = fa[n_front * f_k + n_bottom * f_j + i + offset];
  a_bottom_back = fa[n_back * f_k + n_bottom * f_j + i + offset];

  double tmp = 0.0;
  fxx = (a_right - 2 * a_middle + a_left) / hx / hx;
  fyy = (a_bottom - 2 * a_middle + a_top) / hy / hy;
  fzz = (a_back - 2 * a_middle + a_front) / hz / hz;
  fxy = ((a_right_bottom - a_left_bottom) - (a_right_top - a_left_top)) / 2.0 / hx / 2.0 / hy;
  fxz = ((a_right_back - a_left_back) - (a_right_front - a_left_front)) / 2.0 / hx / 2.0 / hz;
  fyz = ((a_bottom_back - a_top_back) - (a_bottom_front - a_top_front)) / 2.0 / hy / 2.0 / hz;
  // int o = orie[k * nx * ny + j * nx + i];
  // if (o % 12 == n)
  //{
  tmp += lambda[0 * 3 + 0 + n * 3 * 3] * fxx;
  tmp += lambda[0 * 3 + 1 + n * 3 * 3] * fxy;
  tmp += lambda[0 * 3 + 2 + n * 3 * 3] * fxz;
  tmp += lambda[1 * 3 + 0 + n * 3 * 3] * fxy;
  tmp += lambda[1 * 3 + 1 + n * 3 * 3] * fyy;
  tmp += lambda[1 * 3 + 2 + n * 3 * 3] * fyz;
  tmp += lambda[2 * 3 + 0 + n * 3 * 3] * fxz;
  tmp += lambda[2 * 3 + 1 + n * 3 * 3] * fyz;
  tmp += lambda[2 * 3 + 2 + n * 3 * 3] * fzz;
  f[k * nx * ny + j * nx + i] = tmp;
  //}
}

__global__ void
top_bottom_pack(double *field, double *fields_top, double *fields_bottom, double *fieldr_front, double *fieldr_back)
{
  int k_s, k_e;
  int f_j, f_k;
  int tb_j, tb_k;
  int t_ofst, b_ofst, i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  k_s = nghost + 2;
  k_e = nz + (nghost + 2);
  f_j = nx;
  f_k = nx * ny;
  tb_j = nx;
  tb_k = nx * (nghost + 2);
  t_ofst = nx * nghost;
  b_ofst = nx * (ny - nghost - (nghost + 2));

  if (k < k_s)
  {
    if (front >= 0)
    {
      fields_top[tb_k * k + tb_j * j + i] = fieldr_front[f_k * k + f_j * j + i + t_ofst];
      fields_bottom[tb_k * k + tb_j * j + i] = fieldr_front[f_k * k + f_j * j + i + b_ofst];
    }
  }
  else if (k >= k_e)
  {
    if (back >= 0)
    {
      fields_top[tb_k * k + tb_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * j + i + t_ofst];
      fields_bottom[tb_k * k + tb_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * j + i + b_ofst];
    }
  }
  else
  {
    fields_top[tb_k * k + tb_j * j + i] = field[f_k * (k - k_s) + f_j * j + i + t_ofst];
    fields_bottom[tb_k * k + tb_j * j + i] = field[f_k * (k - k_s) + f_j * j + i + b_ofst];
  }
}

__global__ void
left_right_pack(double *field, double *fields_left, double *fields_right, double *fieldr_top, double *fieldr_bottom, double *fieldr_front, double *fieldr_back)
{
  int j_s, j_e, k_s, k_e;
  int f_j, f_k;
  int lr_j, lr_k;
  int tb_j, tb_k;
  int l_ofst, r_ofst, i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

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

  if (j < j_s)
  {
    if (top >= 0)
    {
      fields_left[lr_k * k + lr_j * j + i] = fieldr_top[tb_k * k + tb_j * j + i + l_ofst];
      fields_right[lr_k * k + lr_j * j + i] = fieldr_top[tb_k * k + tb_j * j + i + r_ofst];
    }
  }
  else if (j >= j_e)
  {
    if (bottom >= 0)
    {
      fields_left[lr_k * k + lr_j * j + i] = fieldr_bottom[tb_k * k + tb_j * (j - j_e) + i + l_ofst];
      fields_right[lr_k * k + lr_j * j + i] = fieldr_bottom[tb_k * k + tb_j * (j - j_e) + i + r_ofst];
    }
  }
  else
  {
    if (k < k_s)
    {
      if (front >= 0)
      {
        fields_left[lr_k * k + lr_j * j + i] = fieldr_front[f_k * k + f_j * (j - j_s) + i + l_ofst];
        fields_right[lr_k * k + lr_j * j + i] = fieldr_front[f_k * k + f_j * (j - j_s) + i + r_ofst];
      }
    }
    else if (k >= k_e)
    {
      if (back >= 0)
      {
        fields_left[lr_k * k + lr_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * (j - j_s) + i + l_ofst];
        fields_right[lr_k * k + lr_j * j + i] = fieldr_back[f_k * (k - k_e) + f_j * (j - j_s) + i + r_ofst];
      }
    }
    else
    {
      fields_left[lr_k * k + lr_j * j + i] = field[f_k * (k - k_s) + f_j * (j - j_s) + i + l_ofst];
      fields_right[lr_k * k + lr_j * j + i] = field[f_k * (k - k_s) + f_j * (j - j_s) + i + r_ofst];
    }
  }
}

__global__ void
left_right_unpack(double *field, double *fieldr_left, double *fieldr_right)
{
  int f_j, f_k;
  int lr_k, lr_j;
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f_j = nx;
  f_k = nx * ny;
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);

  if (left >= 0)
  {
    field[f_k * k + f_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i + 2];
  }
  if (right >= 0)
  {
    field[f_k * k + f_j * j + i + (nx - nghost)] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i];
  }
}

__global__ void
top_bottom_unpack(double *field, double *fieldr_top, double *fieldr_bottom)
{

  int f_j, f_k;
  int tb_k, tb_j;
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f_j = nx;
  f_k = nx * ny;
  tb_j = nx;
  tb_k = nx * (nghost + 2);

  if (top >= 0)
  {
    field[f_k * k + f_j * j + i] = fieldr_top[tb_k * (k + (nghost + 2)) + tb_j * (j + 2) + i];
  }
  if (bottom >= 0)
  {
    field[f_k * k + f_j * (j + ny - nghost) + i] = fieldr_bottom[tb_k * (k + (nghost + 2)) + tb_j * j + i];
  }
}

__global__ void
front_back_unpack(double *field, double *fieldr_front, double *fieldr_back)
{

  int f_j, f_k;
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f_j = nx;
  f_k = nx * ny;

  if (front >= 0)
  {
    field[f_k * k + f_j * j + i] = fieldr_front[f_k * (k + 2) + f_j * j + i];
  }
  if (back >= 0)
  {
    field[f_k * (k + nz - nghost) + f_j * j + i] = fieldr_back[f_k * k + f_j * j + i];
  }
}

__global__ void
left_right_enlarge(double *fielde_left, double *fielde_right, double *fieldr_left, double *fieldr_right)
{
  int j_s, j_e, k_s, k_e;
  int elr_k, elr_j;
  int lr_k, lr_j;
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

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
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);

  if (left >= 0)
  {
    if (front >= 0 && k < k_s)
    {
      if (j >= 2 && j < ny + 2)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + nghost) + i];
      }
      if (top >= 0 && j < nghost + 2)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * k + lr_j * j + i];
      }
      if (bottom >= 0 && j >= ny)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + (nghost + 2)) + i];
      }
    }
    if (k >= k_s && k < k_e)
    {
      if (j >= 2 && j < ny + 2)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * (j + nghost) + i];
      }
      if (top >= 0 && j < nghost + 2)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * j + i];
      }
      if (bottom >= 0 && j >= ny)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * (j + (nghost + 2)) + i];
      }
    }
    if (back >= 0 && k >= k_e)
    {
      if (j >= 2 && j < ny + 2)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + nghost) + i];
      }
      if (top >= 0 && j < nghost + 2)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * j + i];
      }
      if (bottom >= 0 && j >= ny)
      {
        fielde_left[elr_k * k + elr_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i];
      }
    }
  }

  if (right >= 0)
  {
    if (front >= 0 && k < k_s)
    {
      if (j >= 2 && j < ny + 2)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + nghost) + i];
      }
      if (top >= 0 && j < nghost + 2)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * k + lr_j * j + i];
      }
      if (bottom >= 0 && j >= ny)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + (nghost + 2)) + i];
      }
    }
    if (k >= k_s && k < k_e)
    {
      if (j >= 2 && j < ny + 2)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * (j + nghost) + i];
      }
      if (top >= 0 && j < nghost + 2)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * j + i];
      }
      if (bottom >= 0 && j >= ny)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * (j + (nghost + 2)) + i];
      }
    }
    if (back >= 0 && k >= k_e)
    {
      if (j >= 2 && j < ny + 2)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + nghost) + i];
      }
      if (top >= 0 && j < nghost + 2)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * j + i];
      }
      if (bottom >= 0 && j >= ny)
      {
        fielde_right[elr_k * k + elr_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + (nghost + 2)) + i];
      }
    }
  }
}

__global__ void
top_bottom_enlarge(double *fielde_top, double *fielde_bottom, double *fieldr_left, double *fieldr_right, double *fieldr_top, double *fieldr_bottom)
{
  int j_s, j_e, k_s, k_e;
  int etb_k, etb_j;
  int lr_k, lr_j, tb_k, tb_j;
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

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

  etb_j = nx + 4;
  etb_k = (nx + 4) * (nghost + 2);
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);
  tb_j = nx;
  tb_k = nx * (nghost + 2);

  if (top >= 0)
  {
    if (front >= 0 && k < k_s)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_top[tb_k * k + tb_j * j + (i - 2)];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * k + lr_j * j + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * k + lr_j * j + (i - nx)];
      }
    }
    if (k >= k_s && k < k_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_top[tb_k * (k + nghost) + tb_j * j + (i - 2)];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * j + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * j + (i - nx)];
      }
    }
    if (back >= 0 && k >= k_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_top[tb_k * (k + (nghost + 2)) + tb_j * j + (i - 2)];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * j + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_top[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * j + (i - nx)];
      }
    }
  }

  if (bottom >= 0)
  {
    if (front >= 0 && k < k_s)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_bottom[tb_k * k + tb_j * j + (i - 2)];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + (ny + (nghost + 2))) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + (ny + (nghost + 2))) + (i - nx)];
      }
    }
    if (k >= k_s && k < k_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_bottom[tb_k * (k + nghost) + tb_j * j + (i - 2)];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + nghost) + lr_j * (j + (ny + (nghost + 2))) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + nghost) + lr_j * (j + (ny + (nghost + 2))) + (i - nx)];
      }
    }
    if (back >= 0 && k >= k_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_bottom[tb_k * (k + (nghost + 2)) + tb_j * j + (i - 2)];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_left[lr_k * (k + (nghost + 2)) + lr_j * (j + (ny + (nghost + 2))) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_bottom[etb_k * k + etb_j * j + i] = fieldr_right[lr_k * (k + (nghost + 2)) + lr_j * (j + (ny + (nghost + 2))) + (i - nx)];
      }
    }
  }
}

__global__ void
front_back_enlarge(double *fielde_front, double *fielde_back, double *fieldr_left, double *fieldr_right, double *fieldr_top, double *fieldr_bottom, double *fieldr_front, double *fieldr_back)
{
  int j_s, j_e, k_s, k_e;
  int efb_k, efb_j;
  int lr_k, lr_j, tb_k, tb_j, fb_k, fb_j;
  int b_ofst;
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

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

  efb_j = (nx + 4) * (nghost + 2);
  efb_k = nx + 4;
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);
  tb_j = nx;
  tb_k = nx * (nghost + 2);
  fb_j = nx;
  fb_k = nx * ny;
  b_ofst = nz + (nghost + 2);

  if (front >= 0)
  {
    if (top >= 0 && j < j_s)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_top[tb_k * k + tb_j * j + i - 2];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * k + lr_j * j + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * k + lr_j * j + i - nx];
      }
    }
    if (j >= j_s && j < j_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_front[fb_k * k + fb_j * (j - nghost) + i - 2];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + nghost) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + nghost) + i - nx];
      }
    }
    if (bottom >= 0 && j >= j_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_bottom[tb_k * k + tb_j * (j - j_e) + i - 2];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * k + lr_j * (j + (nghost + 2)) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_front[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * k + lr_j * (j + (nghost + 2)) + i - nx];
      }
    }
  }
  if (back >= 0)
  {
    if (top >= 0 && j < j_s)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_top[tb_k * (k + b_ofst) + tb_j * j + i - 2];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * (k + b_ofst) + lr_j * j + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * (k + b_ofst) + lr_j * j + i - nx];
      }
    }
    if (j >= j_s && j < j_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_back[fb_k * k + fb_j * (j - nghost) + i - 2];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * (k + b_ofst) + lr_j * (j + nghost) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * (k + b_ofst) + lr_j * (j + nghost) + i - nx];
      }
    }
    if (bottom >= 0 && j >= j_e)
    {
      if (i >= 2 && i < nx + 2)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_bottom[tb_k * (k + b_ofst) + tb_j * (j - j_e) + i - 2];
      }
      if (left >= 0 && i < nghost + 2)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_left[lr_k * (k + b_ofst) + lr_j * (j + (nghost + 2)) + i];
      }
      if (right >= 0 && i >= nx)
      {
        fielde_back[efb_k * k + efb_j * j + i] = fieldr_right[lr_k * (k + b_ofst) + lr_j * (j + (nghost + 2)) + i - nx];
      }
    }
  }
}

__global__ void
left_right_mu(double *Eu_left, double *Eu_right, double *Ee_left, double *Ee_right)
{
  int i, j, k;
  int mu_i, mu_j, mu_k;
  int e_i, e_j, e_k;
  int l_e, l_mu;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y + 2;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z + 2;

  mu_k = ny;
  e_j = nghost + 2;
  e_k = (nghost + 2) * (ny + 4);

  if (left >= 0)
  {

    if (k == 2 && front < 0)
    {
      l_e = k * e_k + j * e_j + 1;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_left[l_mu] = Ee_left[l_e];
    }
    else if (k == nz + 1 && back < 0)
    {
      l_e = k * e_k + j * e_j + 1;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_left[l_mu] = Ee_left[l_e];
    }
    else
    {
      l_e = k * e_k + j * e_j + 1;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_left[l_mu] = Ee_left[l_e];
    }
  }

  if (right >= 0)
  {

    if (k == 2 && front < 0)
    {
      l_e = k * e_k + j * e_j + 2;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_right[l_mu] = Ee_right[l_e];
    }
    else if (k == nz + 1 && back < 0)
    {
      l_e = k * e_k + j * e_j + 2;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_right[l_mu] = Ee_right[l_e];
    }
    else
    {
      l_e = k * e_k + j * e_j + 2;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_right[l_mu] = Ee_right[l_e];
    }
  }
}

__global__ void
top_bottom_mu(double *Eu_top, double *Eu_bottom, double *Ee_top, double *Ee_bottom)
{
  int i, j, k;
  int mu_i, mu_j, mu_k;
  int e_i, e_j, e_k;
  int l_e, l_mu;

  i = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  j = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x + 2;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z + 2;

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  if (top >= 0)
  {

    if (k == 2 && front < 0)
    {
      l_e = k * e_k + e_j + j;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_top[l_mu] = Ee_top[l_e];
    }
    else if (k == nz + 1 && back < 0)
    {
      l_e = k * e_k + e_j + j;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_top[l_mu] = Ee_top[l_e];
    }
    else
    {
      l_e = k * e_k + e_j + j;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_top[l_mu] = Ee_top[l_e];
    }
  }

  if (bottom >= 0)
  {

    if (k == 2 && front < 0)
    {
      l_e = k * e_k + e_j * nghost + j;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_bottom[l_mu] = Ee_bottom[l_e];
    }
    else if (k == nz + 1 && back < 0)
    {
      l_e = k * e_k + e_j * nghost + j;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_bottom[l_mu] = Ee_bottom[l_e];
    }
    else
    {
      l_e = k * e_k + e_j * nghost + j;
      l_mu = (k - 2) * mu_k + j - 2;
      Eu_bottom[l_mu] = Ee_bottom[l_e];
    }
  }
}

__global__ void
front_back_mu(Dtype *Eu_front, Dtype *Eu_back, Dtype *Ee_front, Dtype *Ee_back)
{
  int i, j, k;
  int mu_i, mu_j, mu_k;
  int e_i, e_j, e_k;
  int l_e, l_mu;

  i = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;
  j = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x + 2;
  k = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y + 2;

  mu_k = nx;
  e_j = nx + 4;
  e_k = (nx + 4) * (nghost + 2);

  if (front >= 0)
  {

    l_e = k * e_k + e_j + j;
    l_mu = (k - 2) * mu_k + j - 2;
    Eu_front[l_mu] = Ee_front[l_e];
  }

  if (back >= 0)
  {

    l_e = k * e_k + e_j * nghost + j;
    l_mu = (k - 2) * mu_k + j - 2;
    Eu_back[l_mu] = Ee_back[l_e];
  }
}

__global__ void
unpack_lr(double *fieldr, double *fieldr_left, double *fieldr_right)
{
  int f_j, f_k;
  int lr_k, lr_j;
  int i, j, k;

  f_j = nx + 2 * 2;
  f_k = (nx + 2 * 2) * (ny + 2 * 2);
  lr_j = nghost + 2;
  lr_k = (nghost + 2) * (ny + (nghost + 2) * 2);

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (left >= 0)
  {
    fieldr[f_k * k + f_j * j + i] = fieldr_left[lr_k * (k + 2) + lr_j * (j + 2) + i];
  }
  if (right >= 0)
  {
    fieldr[f_k * k + f_j * j + i + nx] = fieldr_right[lr_k * (k + 2) + lr_j * (j + 2) + i];
  }
}

__global__ void
unpack_tb(double *fieldr, double *fieldr_top, double *fieldr_bottom)
{
  int f_j, f_k;
  int tb_k, tb_j;
  int i, j, k;
  f_j = nx + 2 * 2;
  f_k = (nx + 2 * 2) * (ny + 2 * 2);
  tb_j = nx;
  tb_k = nx * (nghost + 2);
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (top >= 0)
  {
    fieldr[f_k * k + f_j * j + i + 2] = fieldr_top[tb_k * (k + 2) + tb_j * j + i];
  }
  if (bottom >= 0)
  {
    fieldr[f_k * k + f_j * (j + ny) + i + 2] = fieldr_bottom[tb_k * (k + 2) + tb_j * j + i];
  }
}

__global__ void
unpack_fb(double *fieldr, double *fieldr_front, double *fieldr_back)
{
  int f_j, f_k;
  int fb_j, fb_k;
  int i, j, k;

  f_j = nx + 2 * 2;
  f_k = (nx + 2 * 2) * (ny + 2 * 2);
  fb_j = nx;
  fb_k = nx * ny;
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (front >= 0)
  {
    fieldr[f_k * k + f_j * (j + 2) + (i + 2)] = fieldr_front[fb_k * k + fb_j * j + i];
  }

  if (back >= 0)
  {
    fieldr[f_k * (k + nz) + f_j * (j + 2) + (i + 2)] = fieldr_back[fb_k * k + fb_j * j + i];
  }
}

__global__ void
unpack_all(double *fieldr, double *field)
{
  int f_j, f_k;
  int i, j, k;
  f_j = nx + 2 * 2;
  f_k = (nx + 2 * 2) * (ny + 2 * 2);
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  fieldr[f_k * (k + 2) + f_j * (j + 2) + (i + 2)] = field[nx * ny * k + nx * j + i];
}

void transfer(void)
{
  int n;

  int threads_x, threads_y, threads_z;

#ifdef SCLETD_DEBUG
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    hipMemcpy(ac[n].fieldEs_front, &fieldE[n * offset + nx * ny * nghost], sizeof(Dtype) * fb_size, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEs_back, &fieldE[n * offset + nx * ny * (nz - nghost - (nghost + 2))], sizeof(Dtype) * fb_size, hipMemcpyDeviceToHost);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(ch[n].fieldCIs_front, &fieldCI[n * offset + nx * ny * nghost], sizeof(Dtype) * fb_size, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIs_back, &fieldCI[n * offset + nx * ny * (nz - nghost - (nghost + 2))], sizeof(Dtype) * fb_size, hipMemcpyDeviceToHost);
  }
#endif

#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_Memcpy_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    MPI_Startall(4, ac[n].ireq_front_back_fieldE);
    MPI_Waitall(4, ac[n].ireq_front_back_fieldE, status);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    MPI_Startall(4, ch[n].ireq_front_back_fieldCI);
    MPI_Waitall(4, ch[n].ireq_front_back_fieldCI, status);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_MPI_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    hipMemcpy(fieldEr_front + n * fb_size, ac[n].fieldEr_front, sizeof(Dtype) * fb_size, hipMemcpyHostToDevice);
    hipMemcpy(fieldEr_back + n * fb_size, ac[n].fieldEr_back, sizeof(Dtype) * fb_size, hipMemcpyHostToDevice);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(fieldCIr_front + n * fb_size, ch[n].fieldCIr_front, sizeof(Dtype) * fb_size, hipMemcpyHostToDevice);
    hipMemcpy(fieldCIr_back + n * fb_size, ch[n].fieldCIr_back, sizeof(Dtype) * fb_size, hipMemcpyHostToDevice);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_Memcpy_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    dim3 blocks_tb_pack(nx / THREADS_PER_BLOCK_X, (nghost + 2) / THREADS_PER_BLOCK_Y, (nz + (nghost + 2) * 2) / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_pack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_pack, blocks_tb_pack, threads_tb_pack, 0, 0,
                       fieldE + n * offset, fieldEs_top + n * tb_size, fieldEs_bottom + n * tb_size, fieldEr_front + n * fb_size, fieldEr_back + n * fb_size);
  }

#if 1
#if 0
  for (n = 0; n < nch; n++)
  {
    dim3 blocks_tb_pack(nx / THREADS_PER_BLOCK_X, (nghost + 2) / THREADS_PER_BLOCK_Y, (nz + (nghost + 2) * 2) / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_pack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_pack, blocks_tb_pack, threads_tb_pack, 0, 0,
                       fieldCI + n * offset, fieldCIs_top + n * tb_size, fieldCIs_bottom + n * tb_size, fieldCIr_front + n * fb_size, fieldCIr_back + n * fb_size);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_pack_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    hipMemcpy(ac[n].fieldEs_top, fieldEs_top + n * tb_size, sizeof(Dtype) * tb_size, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEs_bottom, fieldEs_bottom + n * tb_size, sizeof(Dtype) * tb_size, hipMemcpyDeviceToHost);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(ch[n].fieldCIs_top, fieldCIs_top + n * tb_size, sizeof(Dtype) * tb_size, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIs_bottom, fieldCIs_bottom + n * tb_size, sizeof(Dtype) * tb_size, hipMemcpyDeviceToHost);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_Memcpy_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    MPI_Startall(4, ac[n].ireq_top_bottom_fieldE);
    MPI_Waitall(4, ac[n].ireq_top_bottom_fieldE, status);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    MPI_Startall(4, ch[n].ireq_top_bottom_fieldCI);
    MPI_Waitall(4, ch[n].ireq_top_bottom_fieldCI, status);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_MPI_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    hipMemcpy(fieldEr_top + n * tb_size, ac[n].fieldEr_top, sizeof(Dtype) * tb_size, hipMemcpyHostToDevice);
    hipMemcpy(fieldEr_bottom + n * tb_size, ac[n].fieldEr_bottom, sizeof(Dtype) * tb_size, hipMemcpyHostToDevice);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(fieldCIr_top + n * tb_size, ch[n].fieldCIr_top, sizeof(Dtype) * tb_size, hipMemcpyHostToDevice);
    hipMemcpy(fieldCIr_bottom + n * tb_size, ch[n].fieldCIr_bottom, sizeof(Dtype) * tb_size, hipMemcpyHostToDevice);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_Memcpy_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    dim3 blocks_lr_pack((nghost + 2) / THREADS_PER_BLOCK_X, (ny + (nghost + 2) * 2) / THREADS_PER_BLOCK_Y, (nz + (nghost + 2) * 2) / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_pack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_pack, blocks_lr_pack, threads_lr_pack, 0, 0,
                       fieldE + n * offset, fieldEs_left + n * lr_size, fieldEs_right + n * lr_size,
                       fieldEr_top + n * tb_size, fieldEr_bottom + n * tb_size, fieldEr_front + n * fb_size,
                       fieldEr_back + n * fb_size);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    dim3 blocks_lr_pack((nghost + 2) / THREADS_PER_BLOCK_X, (ny + (nghost + 2) * 2) / THREADS_PER_BLOCK_Y, (nz + (nghost + 2) * 2) / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_pack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_pack, blocks_lr_pack, threads_lr_pack, 0, 0,
                       fieldCI + n * offset, fieldCIs_left + n * lr_size, fieldCIs_right + n * lr_size,
                       fieldCIr_top + n * tb_size, fieldCIr_bottom + n * tb_size, fieldCIr_front + n * fb_size,
                       fieldCIr_back + n * fb_size);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_pack_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    hipMemcpy(ac[n].fieldEs_left, fieldEs_left + n * lr_size, sizeof(Dtype) * lr_size, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEs_right, fieldEs_right + n * lr_size, sizeof(Dtype) * lr_size, hipMemcpyDeviceToHost);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(ch[n].fieldCIs_left, fieldCIs_left + n * lr_size, sizeof(Dtype) * lr_size, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIs_right, fieldCIs_right + n * lr_size, sizeof(Dtype) * lr_size, hipMemcpyDeviceToHost);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_Memcpy_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    MPI_Startall(4, ac[n].ireq_left_right_fieldE);
    MPI_Waitall(4, ac[n].ireq_left_right_fieldE, status);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    MPI_Startall(4, ch[n].ireq_left_right_fieldCI);
    MPI_Waitall(4, ch[n].ireq_left_right_fieldCI, status);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_MPI_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    hipMemcpy(fieldEr_left + n * lr_size, ac[n].fieldEr_left, sizeof(Dtype) * lr_size, hipMemcpyHostToDevice);
    hipMemcpy(fieldEr_right + n * lr_size, ac[n].fieldEr_right, sizeof(Dtype) * lr_size, hipMemcpyHostToDevice);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(fieldCIr_left + n * lr_size, ch[n].fieldCIr_left, sizeof(Dtype) * lr_size, hipMemcpyHostToDevice);
    hipMemcpy(fieldCIr_right + n * lr_size, ch[n].fieldCIr_right, sizeof(Dtype) * lr_size, hipMemcpyHostToDevice);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_Memcpy_time += timer;
#endif

  for (n = 0; n < nac; n++)
  {

#ifdef SCLETD_DEBUG
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    dim3 blocks_lr_enlarge((nghost + 2) / THREADS_PER_BLOCK_X, (ny + 4) / THREADS_PER_BLOCK_Y, (nz + 4) / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_enlarge(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_enlarge, blocks_lr_enlarge, threads_lr_enlarge, 0, 0,
                       fieldEe_left + n * e_lr, fieldEe_right + n * e_lr,
                       fieldEr_left + n * lr_size, fieldEr_right + n * lr_size);

    dim3 blocks_tb_enlarge((nx + 4) / THREADS_PER_BLOCK_X, (nghost + 2) / THREADS_PER_BLOCK_Y, (nz + 4) / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_enlarge(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_enlarge, blocks_tb_enlarge, threads_tb_enlarge, 0, 0,
                       fieldEe_top + n * e_tb, fieldEe_bottom + n * e_tb,
                       fieldEr_left + n * lr_size, fieldEr_right + n * lr_size,
                       fieldEr_top + n * tb_size, fieldEr_bottom + n * tb_size);

    dim3 blocks_fb_enlarge((nx + 4) / THREADS_PER_BLOCK_X, (ny + 4) / THREADS_PER_BLOCK_Y, (nghost + 2) / THREADS_PER_BLOCK_Z);
    dim3 threads_fb_enlarge(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(front_back_enlarge, blocks_fb_enlarge, threads_fb_enlarge, 0, 0,
                       fieldEe_front + n * e_fb, fieldEe_back + n * e_fb,
                       fieldEr_left + n * lr_size, fieldEr_right + n * lr_size,
                       fieldEr_top + n * tb_size, fieldEr_bottom + n * tb_size, fieldEr_front + n * fb_size,
                       fieldEr_back + n * fb_size);

#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    trans_enlarge_time += timer;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    threads_x = 1;
    dim3 blocks_lr_mu(1 / threads_x, ny / THREADS_PER_BLOCK_Y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_mu(threads_x, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_mu, blocks_lr_mu, threads_lr_mu, 0, 0,
                       fieldEu_left + n * u_lr, fieldEu_right + n * u_lr,
                       fieldEe_left + n * e_lr, fieldEe_right + n * e_lr);

    threads_y = 1;
    dim3 blocks_tb_mu(nx / THREADS_PER_BLOCK_X, 1 / threads_y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_mu(THREADS_PER_BLOCK_X, threads_y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_mu, blocks_tb_mu, threads_tb_mu, 0, 0,
                       fieldEu_top + n * u_tb, fieldEu_bottom + n * u_tb,
                       fieldEe_top + n * e_tb, fieldEe_bottom + n * e_tb);

    threads_z = 1;
    dim3 blocks_fb_mu(nx / THREADS_PER_BLOCK_X, ny / THREADS_PER_BLOCK_Y, 1 / threads_z);
    dim3 threads_fb_mu(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, threads_z);
    hipLaunchKernelGGL(front_back_mu, blocks_fb_mu, threads_fb_mu, 0, 0,
                       fieldEu_front + n * u_fb, fieldEu_back + n * u_fb,
                       fieldEe_front + n * e_fb, fieldEe_back + n * e_fb);

#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    trans_mu_time += timer;
#endif
  }

#ifdef SCLETD_DEBUG
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif
#if 0
  for (n = 0; n < nch; n++)
  {
    dim3 blocks_lr_enlarge((nghost + 2) / THREADS_PER_BLOCK_X, (ny + 4) / THREADS_PER_BLOCK_Y, (nz + 4) / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_enlarge(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_enlarge, blocks_lr_enlarge, threads_lr_enlarge, 0, 0,
                       fieldCIe_left + n * e_lr, fieldCIe_right + n * e_lr,
                       fieldCIr_left + n * lr_size, fieldCIr_right + n * lr_size);

    dim3 blocks_tb_enlarge((nx + 4) / THREADS_PER_BLOCK_X, (nghost + 2) / THREADS_PER_BLOCK_Y, (nz + 4) / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_enlarge(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_enlarge, blocks_tb_enlarge, threads_tb_enlarge, 0, 0,
                       fieldCIe_top + n * e_tb, fieldCIe_bottom + n * e_tb,
                       fieldCIr_left + n * lr_size, fieldCIr_right + n * lr_size,
                       fieldCIr_top + n * tb_size, fieldCIr_bottom + n * tb_size);

    dim3 blocks_fb_enlarge((nx + 4) / THREADS_PER_BLOCK_X, (ny + 4) / THREADS_PER_BLOCK_Y, (nghost + 2) / THREADS_PER_BLOCK_Z);
    dim3 threads_fb_enlarge(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(front_back_enlarge, blocks_fb_enlarge, threads_fb_enlarge, 0, 0,
                       fieldCIe_front + n * e_fb, fieldCIe_back + n * e_fb,
                       fieldCIr_left + n * lr_size, fieldCIr_right + n * lr_size,
                       fieldCIr_top + n * tb_size, fieldCIr_bottom + n * tb_size,
                       fieldCIr_front + n * fb_size, fieldCIr_back + n * fb_size);
  }
#endif
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_enlarge_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif
#if 0
  for (n = 0; n < nac; n++)
  {
    hipMemcpy(ac[n].fieldEe_left, fieldEe_left + n * e_lr, sizeof(Dtype) * e_lr, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEe_right, fieldEe_right + n * e_lr, sizeof(Dtype) * e_lr, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEe_top, fieldEe_top + n * e_tb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEe_bottom, fieldEe_bottom + n * e_tb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEe_front, fieldEe_front + n * e_fb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
    hipMemcpy(ac[n].fieldEe_back, fieldEe_back + n * e_fb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
  }
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(ch[n].fieldCIe_left, fieldCIe_left + n * e_lr, sizeof(Dtype) * e_lr, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIe_right, fieldCIe_right + n * e_lr, sizeof(Dtype) * e_lr, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIe_top, fieldCIe_top + n * e_tb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIe_bottom, fieldCIe_bottom + n * e_tb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIe_front, fieldCIe_front + n * e_fb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
    hipMemcpy(ch[n].fieldCIe_back, fieldCIe_back + n * e_fb, sizeof(Dtype) * e_tb, hipMemcpyDeviceToHost);
  }
  for (n = 0; n < nch; n++)
  {
    ch_mu(n, ch[n].LE, ch[n].fieldCImu_left, ch[n].fieldCImu_right, ch[n].fieldCImu_top, ch[n].fieldCImu_bottom, ch[n].fieldCImu_front, ch[n].fieldCImu_back,
          ch[n].fieldCIe_left, ch[n].fieldCIe_right, ch[n].fieldCIe_top, ch[n].fieldCIe_bottom, ch[n].fieldCIe_front, ch[n].fieldCIe_back);
  }
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(fieldCImu_left + n * u_lr, ch[n].fieldCImu_left, sizeof(Dtype) * u_lr, hipMemcpyHostToDevice);
    hipMemcpy(fieldCImu_right + n * u_lr, ch[n].fieldCImu_right, sizeof(Dtype) * u_lr, hipMemcpyHostToDevice);
    hipMemcpy(fieldCImu_top + n * u_tb, ch[n].fieldCImu_top, sizeof(Dtype) * u_tb, hipMemcpyHostToDevice);
    hipMemcpy(fieldCImu_bottom + n * u_tb, ch[n].fieldCImu_bottom, sizeof(Dtype) * u_tb, hipMemcpyHostToDevice);
    hipMemcpy(fieldCImu_front + n * u_fb, ch[n].fieldCImu_front, sizeof(Dtype) * u_fb, hipMemcpyHostToDevice);
    hipMemcpy(fieldCImu_back + n * u_fb, ch[n].fieldCImu_back, sizeof(Dtype) * u_fb, hipMemcpyHostToDevice);
  }
#endif
/*
    threads_x = 1;
    dim3 blocks_lr_mu (1 / threads_x, ny / THREADS_PER_BLOCK_Y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_mu (threads_x, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL (ch_left_right_mu, blocks_lr_mu, threads_lr_mu, 0, 0,\
                        fieldCImu_left + n * u_lr, fieldCImu_right + n * u_lr, fieldCIe_left, fieldCIe_right);

    threads_y = 1;
    dim3 blocks_tb_mu (nx / THREADS_PER_BLOCK_X, 1 / threads_y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_mu (THREADS_PER_BLOCK_X, threads_y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL (ch_top_bottom_mu, blocks_tb_mu, threads_tb_mu, 0, 0,\
                        fieldCImu_top + n * u_tb, fieldCImu_bottom + n * u_tb, fieldCIe_top, fieldCIe_bottom);

    threads_z = 1;
    dim3 blocks_fb_mu (nx / THREADS_PER_BLOCK_X, ny / THREADS_PER_BLOCK_Y, 1 / threads_z);
    dim3 threads_fb_mu (THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, threads_z);
    hipLaunchKernelGGL (ch_front_back_mu, blocks_fb_mu, threads_fb_mu, 0, 0,\
                        fieldCImu_front + n * u_fb, fieldCImu_back + n * u_fb, fieldCIe_front, fieldCIe_back);
*/
#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_mu_time += timer;
  hipEventRecord(st, NULL);
  hipEventSynchronize(st);
#endif

  for (n = 0; n < nac; n++)
  {
    threads_x = nghost;
    dim3 blocks_lr_unpack(nghost / threads_x, ny / THREADS_PER_BLOCK_Y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_unpack(threads_x, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_unpack, blocks_lr_unpack, threads_lr_unpack, 0, 0,
                       fieldE + n * offset, fieldEr_left + n * lr_size, fieldEr_right + n * lr_size);
    threads_y = nghost;
    dim3 blocks_tb_unpack(nx / THREADS_PER_BLOCK_X, nghost / threads_y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_unpack(THREADS_PER_BLOCK_X, threads_y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_unpack, blocks_tb_unpack, threads_tb_unpack, 0, 0,
                       fieldE + n * offset, fieldEr_top + n * tb_size, fieldEr_bottom + n * tb_size);

    threads_z = nghost;
    dim3 blocks_fb_unpack(nx / THREADS_PER_BLOCK_X, ny / THREADS_PER_BLOCK_Y, nghost / threads_z);
    dim3 threads_fb_unpack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, threads_z);
    hipLaunchKernelGGL(front_back_unpack, blocks_fb_unpack, threads_fb_unpack, 0, 0,
                       fieldE + n * offset, fieldEr_front + n * fb_size, fieldEr_back + n * fb_size);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    threads_x = nghost;
    dim3 blocks_lr_unpack(nghost / threads_x, ny / THREADS_PER_BLOCK_Y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_lr_unpack(threads_x, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(left_right_unpack, blocks_lr_unpack, threads_lr_unpack, 0, 0,
                       fieldCI + n * offset, fieldCIr_left + n * lr_size, fieldCIr_right + n * lr_size);
    threads_y = nghost;
    dim3 blocks_tb_unpack(nx / THREADS_PER_BLOCK_X, nghost / threads_y, nz / THREADS_PER_BLOCK_Z);
    dim3 threads_tb_unpack(THREADS_PER_BLOCK_X, threads_y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(top_bottom_unpack, blocks_tb_unpack, threads_tb_unpack, 0, 0,
                       fieldCI + n * offset, fieldCIr_top + n * tb_size, fieldCIr_bottom + n * tb_size);

    threads_z = nghost;
    dim3 blocks_fb_unpack(nx / THREADS_PER_BLOCK_X, ny / THREADS_PER_BLOCK_Y, nghost / threads_z);
    dim3 threads_fb_unpack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, threads_z);
    hipLaunchKernelGGL(front_back_unpack, blocks_fb_unpack, threads_fb_unpack, 0, 0,
                       fieldCI + n * offset, fieldCIr_front + n * fb_size, fieldCIr_back + n * fb_size);
  }
#endif
  for (n = 0; n < nac; n++)
  {
    threads_x = nghost;
    dim3 blocks3((nghost + 2) / threads_x, (ny + 2 * 2) / THREADS_PER_BLOCK_Y, (nz + 2 * 2) / THREADS_PER_BLOCK_Z);
    dim3 threads3(threads_x, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(unpack_lr, blocks3, threads3, 0, 0,
                       fieldEr + n * offset_Er, fieldEr_left, fieldEr_right);

    threads_y = nghost;
    dim3 blocks4((nx) / THREADS_PER_BLOCK_X, (nghost + 2) / threads_y, (nz + 2 * 2) / THREADS_PER_BLOCK_Z);
    dim3 threads4(THREADS_PER_BLOCK_X, threads_y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(unpack_tb, blocks4, threads4, 0, 0,
                       fieldEr + n * offset_Er, fieldEr_top, fieldEr_bottom);

    threads_z = nghost;
    dim3 blocks5((nx) / THREADS_PER_BLOCK_X, (ny) / THREADS_PER_BLOCK_Y, (nghost + 2) / threads_z);
    dim3 threads5(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, threads_z);
    hipLaunchKernelGGL(unpack_fb, blocks5, threads5, 0, 0,
                       fieldEr + n * offset_Er, fieldEr_front, fieldEr_back);

    dim3 blocks6((nx) / THREADS_PER_BLOCK_X, (ny) / THREADS_PER_BLOCK_Y, (nz) / THREADS_PER_BLOCK_Z);
    dim3 threads6(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
    hipLaunchKernelGGL(unpack_all, blocks6, threads6, 0, 0,
                       fieldEr + n * offset_Er, fieldE + n * offset);
  }

#ifdef SCLETD_DEBUG
  hipEventRecord(ed, NULL);
  hipEventSynchronize(ed);
  hipEventElapsedTime(&timer, st, ed);
  trans_unpack_time += timer;
#endif
#endif
}

__global__ void
ac_calc_F1_dev(Dtype *f1, Dtype *fieldE, Dtype *fieldCI, Dtype *fieldEu_left, Dtype *fieldEu_right,
               Dtype *fieldEu_top, Dtype *fieldEu_bottom, Dtype *fieldEu_front, Dtype *fieldEu_back,
               Dtype LE, Dtype KE, Dtype epn2, int n)
{
  int i, j, k;
  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double f_left, f_right, f_top, f_bottom, f_front, f_back;
  double c0_bcc, c1_bcc, c2_bcc;
  double c0_hcp, c1_hcp, c2_hcp;
  double u0, u1, u2, u3, un;
  double u4, u5, u6, u7, u8, u9, u10, u11, u12, u13;
  double GBCC, GHCP, dpfeta1, dpfeta2, dqfeta1, dqfeta2, dfeta1, dpfeta11, dpfeta22, dqfeta11, dqfeta22, dfeta2; // G

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  u0 = fieldE[k * nx * ny + j * nx + i];
  u1 = fieldE[k * nx * ny + j * nx + i + offset];
  u2 = 1.0 - u0;
  u3 = 1.0 - u1;
  u4 = fieldE[k * nx * ny + j * nx + i + 2 * offset];
  u5 = fieldE[k * nx * ny + j * nx + i + 3 * offset];
  u6 = fieldE[k * nx * ny + j * nx + i + 4 * offset];
  u7 = fieldE[k * nx * ny + j * nx + i + 5 * offset];

 c0_hcp = 0.1057;
 c1_hcp = 0.0223;
 // c0_hcp = 0.102;
 // c1_hcp = 0.036;
  c2_hcp = 1.0 - c0_hcp - c1_hcp;

 c0_bcc = 0.102;
 c1_bcc = 0.036;
 // c0_bcc = 0.1057;
 // c1_bcc = 0.0223;
  c2_bcc = 1.0 - c0_bcc - c1_bcc;

  //if (c0 < 0.0 || c1 < 0.0)
  //{
    // printf("bad!\n");
  //}

  GBCC = c2_bcc * GBCCTI + c0_bcc * GBCCAL + c1_bcc * GHSERV + R_a * T_a * (c2_bcc * log(c2_bcc) + c0_bcc * log(c0_bcc) + c1_bcc * log(c1_bcc)) + c0_bcc * c1_bcc * (BL12_0 + BL12_1 * (c0_bcc - c1_bcc)) + c2_bcc * c1_bcc * (BL32_0 + BL32_1 * (c2_bcc - c1_bcc) + BL32_2 * (c2_bcc - c1_bcc) * (c2_bcc - c1_bcc)) + c0_bcc * c2_bcc * (BL13_0 + BL13_1 * (c0_bcc - c2_bcc) + BL13_2 * (c0_bcc - c2_bcc) * (c0_bcc - c2_bcc)) + c0_bcc * c1_bcc * c2_bcc * BL132_0;

  GHCP = c2_hcp * GHSERTI + c0_hcp * GHCPAL + c1_hcp * GHCPV + R_a * T_a * (c2_hcp * log(c2_hcp) + c0_hcp * log(c0_hcp) + c1_hcp * log(c1_hcp)) + c0_hcp * c1_hcp * (HL12_0 + HL12_1 * (c0_hcp - c1_hcp)) + c2_hcp * c1_hcp * HL32_0 + c0_hcp * c2_hcp * (HL13_0 + HL13_1 * (c0_hcp - c2_hcp) + HL13_2 * (c0_hcp - c2_hcp) * (c0_hcp - c2_hcp)) + c0_hcp * c1_hcp * c2_hcp * (c0_hcp * HL132_0 + c2_hcp * HL132_1 + c1_hcp * HL132_2);
  /*
    dpfeta1 = 30.0*u2*u2*(1-u2)*(1-u2);
          dqfeta1 = 1.0*wmega*u2*(1-u2)*(1-2.0*u2);
    dpfeta2 = 30.0*u0*u0*(1-u0)*(1-u0);
          dqfeta2 = 1.0*wmega*u0*(1-u0)*(1-2.0*u0);

    dpfeta11 = 30.0*u3*u3*(1-u3)*(1-u3);
          dqfeta11 = 1.0*wmega*u3*(1-u3)*(1-2.0*u3);
    dpfeta22 = 30.0*u1*u1*(1-u1)*(1-u1);
          dqfeta22 = 1.0*wmega*u1*(1-u1)*(1-2.0*u1);

    dfeta1 = -(GBCC*dpfeta1+dqfeta1)+(GHCP*dpfeta2+dqfeta2);
    dfeta2 = -(GBCC*dpfeta11+dqfeta11)+(GHCP*dpfeta22+dqfeta22);

    if (n == 0)
    {
      f1[k * nx * ny + j * nx + i] = LE * KE * u0;
      f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u0 * u1 * u1) / gnormal;
    }
    else
    {
      f1[k * nx * ny + j * nx + i] = LE * KE * u1;
      f1[k * nx * ny + j * nx + i] -= LE * (dfeta2 + 4.0 * ALPHA * u1 * u0 * u0) / gnormal;
    }

  */
  if (n == 0)
  {
    dpfeta1 = 30.0 * u0 * u0 * (1 - u0) * (1 - u0);
    //dqfeta1 = 2.0 * wmega * u0 * (1 - u0) * (1 - 2.0 * u0);
    dqfeta1 = 4.0 * wmega * u0 * (1 - u0) * (2 - u0);
    dfeta1 = -(GBCC * dpfeta1) + (GHCP * dpfeta1) + dqfeta1;
    //dfeta1 = -(GHCP * dpfeta1) + (GBCC * dpfeta1) + dqfeta1;
    f1[k * nx * ny + j * nx + i] = LE * KE * u0;
    // f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u0 * u1 * u1) / gnormal;
    f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 / gnormal + 4.0 * ALPHA * (u0 * u1 * u1 + u0 * u4 * u4 + u0 * u5 * u5 + u0 * u6 * u6 + u0 * u7 * u7));
  }
  else if (n == 1)
  {
    dpfeta1 = 30.0 * u1 * u1 * (1 - u1) * (1 - u1);
    //dqfeta1 = 2.0 * wmega * u1 * (1 - u1) * (1 - 2.0 * u1);
    dqfeta1 = 4.0 * wmega * u1 * (1 - u1) * (2 - u1);
    dfeta1 = -(GBCC * dpfeta1) + (GHCP * dpfeta1) + dqfeta1;
    //dfeta1 = -(GHCP * dpfeta1) + (GBCC * dpfeta1) + dqfeta1;
    f1[k * nx * ny + j * nx + i] = LE * KE * u1;
    // f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u1 * u0 * u0) / gnormal;
    f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 / gnormal + 4.0 * ALPHA * (u1 * u0 * u0 + u1 * u4 * u4 + u1 * u5 * u5 + u1 * u6 * u6 + u1 * u7 * u7));
  }
  else if (n == 2)
  {
    dpfeta1 = 30.0 * u4 * u4 * (1 - u4) * (1 - u4);
    dqfeta1 = 4.0 * wmega * u4 * (1 - u4) * (2 - u4);
    //dqfeta1 = 2.0 * wmega * u4 * (1 - u4) * (1 - 2.0 * u4);
    dfeta1 = -(GBCC * dpfeta1) + (GHCP * dpfeta1) + dqfeta1;
    //dfeta1 = -(GHCP * dpfeta1) + (GBCC * dpfeta1) + dqfeta1;
    f1[k * nx * ny + j * nx + i] = LE * KE * u4;
    // f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u1 * u0 * u0) / gnormal;
    f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 / gnormal + 4.0 * ALPHA * (u4 * u0 * u0 + u4 * u1 * u1 + u4 * u5 * u5 + u4 * u6 * u6 + u4 * u7 * u7));
  }
  else if (n == 3)
  {
    dpfeta1 = 30.0 * u5 * u5 * (1 - u5) * (1 - u5);
    dqfeta1 = 4.0 * wmega * u5 * (1 - u5) * (2 - u5);
    //dqfeta1 = 2.0 * wmega * u5 * (1 - u5) * (1 - 2.0 * u5);
    dfeta1 = -(GBCC * dpfeta1) + (GHCP * dpfeta1) + dqfeta1;
    //dfeta1 = -(GHCP * dpfeta1) + (GBCC * dpfeta1) + dqfeta1;
    f1[k * nx * ny + j * nx + i] = LE * KE * u5;
    // f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u1 * u0 * u0) / gnormal;
    f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 / gnormal + 4.0 * ALPHA * (u5 * u0 * u0 + u5 * u1 * u1 + u5 * u4 * u4 + u5 * u6 * u6 + u5 * u7 * u7));
  }
  else if (n == 4)
  {
    dpfeta1 = 30.0 * u6 * u6 * (1 - u6) * (1 - u6);
    dqfeta1 = 4.0 * wmega * u6 * (1 - u6) * (2 - u6);
    //dqfeta1 = 2.0 * wmega * u5 * (1 - u5) * (1 - 2.0 * u5);
    dfeta1 = -(GBCC * dpfeta1) + (GHCP * dpfeta1) + dqfeta1;
    //dfeta1 = -(GHCP * dpfeta1) + (GBCC * dpfeta1) + dqfeta1;
    f1[k * nx * ny + j * nx + i] = LE * KE * u6;
    // f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u1 * u0 * u0) / gnormal;
    f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 / gnormal + 4.0 * ALPHA * (u6 * u0 * u0 + u6 * u1 * u1 + u6 * u4 * u4 + u6 * u5 * u5 + u6 * u7 * u7));
  }
  else if (n == 5)
  {
    dpfeta1 = 30.0 * u7 * u7 * (1 - u7) * (1 - u7);
    dqfeta1 = 4.0 * wmega * u7 * (1 - u7) * (2 - u7);
    //dqfeta1 = 2.0 * wmega * u5 * (1 - u5) * (1 - 2.0 * u5);
    dfeta1 = -(GBCC * dpfeta1) + (GHCP * dpfeta1) + dqfeta1;
    //dfeta1 = -(GHCP * dpfeta1) + (GBCC * dpfeta1) + dqfeta1;
    f1[k * nx * ny + j * nx + i] = LE * KE * u7;
    // f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 + 4.0 * ALPHA * u1 * u0 * u0) / gnormal;
    f1[k * nx * ny + j * nx + i] -= LE * (dfeta1 / gnormal + 4.0 * ALPHA * (u7 * u0 * u0 + u7 * u1 * u1 + u7 * u4 * u4 + u7 * u6 * u6 + u7 * u5 * u5));
  }

  if (left >= 0)
  {
    if (i == ix1)
    {
      f1[k * nx * ny + j * nx + i] += LE * epn2 * fieldEu_left[k * ny + j] / hx / hx;
    }
  }

  if (right >= 0)
  {
    if (i == ix4 - 1)
    {
      f1[k * nx * ny + j * nx + i] += LE * epn2 * fieldEu_right[k * ny + j] / hx / hx;
    }
  }

  if (top >= 0)
  {
    if (j == iy1)
    {
      f1[k * nx * ny + j * nx + i] += LE * epn2 * fieldEu_top[k * nx + i] / hy / hy;
    }
  }
  if (bottom >= 0)
  {
    if (j == iy4 - 1)
    {
      f1[k * nx * ny + j * nx + i] += LE * epn2 * fieldEu_bottom[k * nx + i] / hy / hy;
    }
  }

  if (front >= 0)
  {
    if (k == iz1)
    {
      f1[k * nx * ny + j * nx + i] += LE * epn2 * fieldEu_front[j * nx + i] / hz / hz;
    }
  }

  if (back >= 0)
  {
    if (k == iz4 - 1)
    {
      f1[k * nx * ny + j * nx + i] += LE * epn2 * fieldEu_back[j * nx + i] / hz / hz;
    }
  }
}

__global__ void
ac_calc_F2_dev(Dtype *f1, Dtype *Elas, Dtype *f2, Dtype LE, Dtype epn2, int *f, Dtype *fieldE, int n)
{
  int i, j, k, l;
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (ELASTIC == 1)
  {
    f1[k * nx * ny + j * nx + i] -= LE * ElasticScale * Elas[k * ny * nx + j * nx + i];
  }
  if (ANISOTROPIC == 1)
  {
    f1[k * nx * ny + j * nx + i] -= epn2 * LE * f2[k * nx * ny + j * nx + i];
  }
  l = k * nx * ny + j * nx + i;
  if (f[l] == 0)
  {
    f1[k * nx * ny + j * nx + i] -= 0.5 * fieldE[k * nx * ny + j * nx + i];
  }
}

__global__ void
update_M_C(int n, double *C, double *M, double *fieldE, double *fieldCI)
{
  int i, j, k;
  double DG1_B, B1_B, DG2_B, B2_B, DG3_B, B3_B,
      DG1_H, B1_H, DG2_H, B2_H, DG3_H, B3_H,
      B1, B2, B3, AA, EE;
  double c0, c1, c2;
  double u0, u1, un;
  double u2, u3, u4, u5, u6, u7, u8, u9, u10, u11;
  C[0] = -1.0e30;
  C[1] = -1.0e30;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  c0 = fieldCI[k * nx * ny + j * nx + i];
  c1 = fieldCI[k * nx * ny + j * nx + i + offset];
  c2 = 1.0 - c0 - c1;
  u0 = fieldE[k * nx * ny + j * nx + i];
  u1 = fieldE[k * nx * ny + j * nx + i + offset];
  u2 = fieldE[k * nx * ny + j * nx + i + 2 * offset];
  u3 = fieldE[k * nx * ny + j * nx + i + 3 * offset];
  un = u0 + u1 + u2 + u3;

  if (un > 1.0)
    un = 1.0;

  // mobility of AL (B1)
  DG1_B = c0 * G1_1 + c1 * G1_2 + c2 * G1_3 + c0 * c1 * G1_12_0 + c2 * c1 * G1_32_0 + c0 * c2 * (G1_13_0 + G1_13_1 * (c0 - c2));
  B1_B = exp(DG1_B / R_a / T_a) / R_a / T_a / Bnormal;

  // mobility of V (B2)
  DG2_B = c0 * G2_1 + c1 * G2_2 + c2 * G2_3 + c0 * c1 * G2_12_0 + c2 * c1 * (G2_32_0 + G2_32_1 * (c2 - c1));
  B2_B = exp(DG2_B / R_a / T_a) / R_a / T_a / Bnormal;

  // mobility of TI (B3)
  DG3_B = c0 * G3_1 + c1 * G3_2 + c2 * G3_3 + c2 * c1 * (G3_32_0 + G3_32_1 * (c2 - c1)) + c0 * c2 * (G3_13_0 + G3_13_1 * (c0 - c2));
  B3_B = exp(DG3_B / R_a / T_a) / R_a / T_a / Bnormal;

  // HCP
  DG1_H = HG1_1;
  B1_H = exp(DG1_H / R_a / T_a) / R_a / T_a / Bnormal;

  // mobility of V (B2)
  DG2_H = HG2_1;
  B2_H = exp(DG2_H / R_a / T_a) / R_a / T_a / Bnormal;

  // mobility of TI (B3)
  DG3_H = HG3_1;
  B3_H = exp(DG3_H / R_a / T_a) / R_a / T_a / Bnormal;

  B1 = exp((DG1_H + un * (DG1_B - DG1_H)) / R_a / T_a) / R_a / T_a / Bnormal;
  B2 = exp((DG2_H + un * (DG2_B - DG2_H)) / R_a / T_a) / R_a / T_a / Bnormal;
  B3 = exp((DG3_H + un * (DG3_B - DG3_H)) / R_a / T_a) / R_a / T_a / Bnormal;

  B1 = -B1 + B1_H + B1_B;
  B2 = -B2 + B2_H + B2_B;
  B3 = -B3 + B3_H + B3_B;

  AA = -(1 - c0) * c0 * B1 * (1 - B3 / B1) + c0 * c1 * B2 * (1 - B3 / B2);
  EE = -(1 - c1) * c1 * B2 * (1 - B3 / B2) + c0 * c1 * B1 * (1 - B3 / B1);

  M[k * nx * ny + j * nx + i] = ((1 - c0) * c0 * B1 + AA * c0);
  M[k * nx * ny + j * nx + i + 1 * offset] = (AA * c1 - c0 * c1 * B2);
  M[k * nx * ny + j * nx + i + 2 * offset] = (EE * c0 - c0 * c1 * B1);
  M[k * nx * ny + j * nx + i + 3 * offset] = ((1 - c1) * c1 * B2 + EE * c1);

  C[0] = MAX(C[0], M[k * ny * nx + j * nx + i + offset * 0]);
  C[0] = MAX(C[0], M[k * ny * nx + j * nx + i + offset * 1]);
  C[1] = MAX(C[1], M[k * ny * nx + j * nx + i + offset * 2]);
  C[1] = MAX(C[1], M[k * ny * nx + j * nx + i + offset * 3]);
}

__global__ void
ch_get_df(int n, double *df1, double *df2, double *fieldE, double *fieldCI)
{
  int i, j, k;
  double c0, c1, c2;
  double u0, u1, u2, u3, un;
  double u4, u5, u6, u7, u8, u9, u10, u11, u12, u13;
  double pfeta1, qfeta1, pfeta2, qfeta2, dGBCC_1, dGBCC_2, dGHCP_1, dGHCP_2, dGall_1, dGall_2;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  c0 = fieldCI[k * nx * ny + j * nx + i];
  c1 = fieldCI[k * nx * ny + j * nx + i + offset];
  c2 = 1 - c0 - c1;

  u0 = fieldE[k * nx * ny + j * nx + i];
  u1 = fieldE[k * nx * ny + j * nx + i + offset];
  u2 = 1.0 - u0;
  u3 = 1.0 - u1;
  u4 = fieldE[k * nx * ny + j * nx + i + 2 * offset];
  u5 = fieldE[k * nx * ny + j * nx + i + 3 * offset];

  if (c0 < 0.0 || c1 < 0.0)
  {
    // printf("bad!\n");
  }

  dGBCC_1 = -GBCCTI + GBCCAL + R_a * T_a * log(c0 / c2) + c1 * (BL12_0 + BL12_1 * (c0 - c1)) + c0 * c1 * BL12_1 - c1 * (BL32_0 + BL32_1 * (c2 - c1) + BL32_2 * (c2 - c1) * (c2 - c1)) - c2 * c1 * (BL32_1 + 2.0 * BL32_2 * (c2 - c1)) + (c2 - c0) * (BL13_0 + BL13_1 * (c0 - c2) + BL13_2 * (c0 - c2) * (c0 - c2)) + c0 * c2 * (2.0 * BL13_1 + 4 * BL13_2 * (c0 - c2)) + BL132_0 * c1 * (c2 - c0);

  dGBCC_2 = -GBCCTI + GHSERV + R_a * T_a * log(c1 / c2) + c0 * (BL12_0 + BL12_1 * (c0 - c1)) - c0 * c1 * BL12_1 + (c2 - c1) * (BL32_0 + BL32_1 * (c2 - c1) + BL32_2 * (c2 - c1) * (c2 - c1)) + c2 * c1 * (-2.0 * BL32_1 - 4 * BL32_2 * (c2 - c1)) - c0 * (BL13_0 + BL13_1 * (c0 - c2) + BL13_2 * (c0 - c2) * (c0 - c2)) + c0 * c2 * (BL13_1 + 2.0 * BL13_2 * (c0 - c2)) + BL132_0 * c0 * (c2 - c1);

  dGHCP_1 = -GHSERTI + GHCPAL + R_a * T_a * log(c0 / c2) + c1 * (HL12_0 + HL12_1 * (c0 - c1)) + c0 * c1 * HL12_1 - c1 * HL32_0 + (c2 - c0) * (HL13_0 + HL13_1 * (c0 - c2) + HL13_2 * (c0 - c2) * (c0 - c2)) + c0 * c2 * (2.0 * HL13_1 + 4 * HL13_2 * (c0 - c2)) + c1 * (c2 - c0) * (c0 * HL132_0 + c2 * HL132_1 + c1 * HL132_2) + c0 * c1 * c2 * (HL132_0 - HL132_1);

  dGHCP_2 = -GHSERTI + GHCPV + R_a * T_a * log(c1 / c2) + c0 * (HL12_0 + HL12_1 * (c0 - c1)) - c0 * c1 * HL12_1 + (c2 - c1) * HL32_0 - c0 * (HL13_0 + HL13_1 * (c0 - c2) + HL13_2 * (c0 - c2) * (c0 - c2)) + c0 * c2 * (HL13_1 + 2.0 * HL13_2 * (c0 - c2)) + c0 * (c2 - c1) * (c0 * HL132_0 + c2 * HL132_1 + c1 * HL132_2) + c0 * c1 * c2 * (HL132_2 - HL132_1);

  //	pfeta1=u2*u2*u2*(10.0-15.0*u2+6.0*u2*u2)+u3*u3*u3*(10.0-15.0*u3+6.0*u3*u3);
  //	pfeta2=u0*u0*u0*(10.0-15.0*u0+6.0*u0*u0)+u1*u1*u1*(10.0-15.0*u1+6.0*u1*u1);

  pfeta2 = u0 * u0 * u0 * (10.0 - 15.0 * u0 + 6.0 * u0 * u0) + u1 * u1 * u1 * (10.0 - 15.0 * u1 + 6.0 * u1 * u1) + u4 * u4 * u4 * (10.0 - 15.0 * u4 + 6.0 * u4 * u4) + u5 * u5 * u5 * (10.0 - 15.0 * u5 + 6.0 * u5 * u5);
  pfeta1 = 1.0 - pfeta2;
  df1[k * nx * ny + j * nx + i] = (pfeta1 * dGBCC_1 + pfeta2 * dGHCP_1) / gnormal;
  df2[k * nx * ny + j * nx + i] = (pfeta1 * dGBCC_2 + pfeta2 * dGHCP_2) / gnormal;
}

__global__ void
ch_calc_F1_A(int n, double *f, double *df1, double *df2, double *fieldCI, double KE)
{
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (n == 0)
  {
    f[k * nx * ny + j * nx + i] = df1[k * nx * ny + j * nx + i] - KE * fieldCI[k * nx * ny + j * nx + i];
  }
  else
  {
    f[k * nx * ny + j * nx + i] = df2[k * nx * ny + j * nx + i] - KE * fieldCI[k * nx * ny + j * nx + i + 1 * offset];
  }
}

__global__ void
ch_calc_F1_B(int n, double *f1, double *f, double *M, double *C)
{
  int i, j, k;
  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double a_left, a_right, a_top, a_bottom, a_front, a_back;
  double a_middle;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

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
    a_front = f[n_front * nx * ny + j * nx + i];
  }
  else
  {
    a_front = 0;
  }
  if (n_back > -1)
  {
    a_back = f[n_back * nx * ny + j * nx + i];
  }
  else
  {
    a_back = 0;
  }
  if (n_top > -1)
  {
    a_top = f[k * nx * ny + n_top * nx + i];
  }
  else
  {
    a_top = 0;
  }
  if (n_bottom > -1)
  {
    a_bottom = f[k * nx * ny + n_bottom * nx + i];
  }
  else
  {
    a_bottom = 0;
  }
  if (n_right > -1)
  {
    a_right = f[k * nx * ny + j * nx + n_right];
  }
  else
  {
    a_right = 0;
  }
  if (n_left > -1)
  {
    a_left = f[k * nx * ny + j * nx + n_left];
  }
  else
  {
    a_left = 0;
  }
  a_middle = f[k * nx * ny + j * nx + i];
  double tmp = 0.0;
  f1[k * nx * ny + j * nx + i] = (a_left + a_right - 2.0 * f[k * nx * ny + j * nx + i]) / hx / hx;
  f1[k * nx * ny + j * nx + i] += (a_top + a_bottom - 2.0 * f[k * nx * ny + j * nx + i]) / hy / hy;
  f1[k * nx * ny + j * nx + i] += (a_front + a_back - 2.0 * f[k * nx * ny + j * nx + i]) / hz / hz;
  f1[k * nx * ny + j * nx + i] = C[n] * f1[k * nx * ny + j * nx + i];

  M[k * nx * ny + j * nx + i + offset * 0] = M[k * nx * ny + j * nx + i + offset * 0] - C[0];
  M[k * nx * ny + j * nx + i + offset * 3] = M[k * nx * ny + j * nx + i + offset * 3] - C[1];
}

__global__ void
ch_calc_F2(int n, double *f2, double *M, double *df)
{
  int i, j, k;
  int n_left, n_right, n_top, n_bottom, n_front, n_back;
  double a_left, a_right, a_top, a_bottom, a_front, a_back, a_middle;
  double m_left, m_right, m_top, m_bottom, m_front, m_back, m_middle;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

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
    a_front = df[n_front * nx * ny + j * nx + i];
  }
  else
  {
    a_front = 0;
  }
  if (n_back > -1)
  {
    a_back = df[n_back * nx * ny + j * nx + i];
  }
  else
  {
    a_back = 0;
  }
  if (n_top > -1)
  {
    a_top = df[k * nx * ny + n_top * nx + i];
  }
  else
  {
    a_top = 0;
  }
  if (n_bottom > -1)
  {
    a_bottom = df[k * nx * ny + n_bottom * nx + i];
  }
  else
  {
    a_bottom = 0;
  }
  if (n_right > -1)
  {
    a_right = df[k * nx * ny + j * nx + n_right];
  }
  else
  {
    a_right = 0;
  }
  if (n_left > -1)
  {
    a_left = df[k * nx * ny + j * nx + n_left];
  }
  else
  {
    a_left = 0;
  }

  if (n_front > -1)
  {
    m_front = M[n_front * nx * ny + j * nx + i];
  }
  else
  {
    m_front = 0;
  }
  if (n_back > -1)
  {
    m_back = M[n_back * nx * ny + j * nx + i];
  }
  else
  {
    m_back = 0;
  }
  if (n_top > -1)
  {
    m_top = M[k * nx * ny + n_top * nx + i];
  }
  else
  {
    m_top = 0;
  }
  if (n_bottom > -1)
  {
    m_bottom = M[k * nx * ny + n_bottom * nx + i];
  }
  else
  {
    m_bottom = 0;
  }
  if (n_right > -1)
  {
    m_right = M[k * nx * ny + j * nx + n_right];
  }
  else
  {
    m_right = 0;
  }
  if (n_left > -1)
  {
    m_left = M[k * nx * ny + j * nx + n_left];
  }
  else
  {
    m_left = 0;
  }
  a_middle = df[k * nx * ny + j * nx + i];
  m_middle = M[k * nx * ny + j * nx + i];

  f2[k * nx * ny + j * nx + i] = (m_left + m_middle) * (a_left - a_middle) / 2.0 / hx / hx;
  f2[k * nx * ny + j * nx + i] += (m_right + m_middle) * (a_right - a_middle) / 2.0 / hx / hx;
  f2[k * nx * ny + j * nx + i] += (m_top + m_middle) * (a_top - a_middle) / 2.0 / hy / hy;
  f2[k * nx * ny + j * nx + i] += (m_bottom + m_middle) * (a_bottom - a_middle) / 2.0 / hy / hy;
  f2[k * nx * ny + j * nx + i] += (m_front + m_middle) * (a_front - a_middle) / 2.0 / hz / hz;
  f2[k * nx * ny + j * nx + i] += (m_back + m_middle) * (a_back - a_middle) / 2.0 / hz / hz;
}

__global__ void
ac_add_F1_F2_dev(int n, Dtype *f1, Dtype *f2, Dtype *f3, Dtype *f4, Dtype *fieldCImu_left, Dtype *fieldCImu_right, Dtype *fieldCImu_top, Dtype *fieldCImu_bottom, Dtype *fieldCImu_front, Dtype *fieldCImu_back, Dtype *C)
{
  int i, j, k;
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f1[k * nx * ny + j * nx + i] = f2[k * nx * ny + j * nx + i];
  f1[k * nx * ny + j * nx + i] += f3[k * nx * ny + j * nx + i];
  f1[k * nx * ny + j * nx + i] += f4[k * nx * ny + j * nx + i];

  if (left >= 0)
  {
    if (i == ix1)
    {
      f1[k * nx * ny + j * nx + i] -= C[n] * fieldCImu_left[k * ny + j] / hx / hx;
    }
  }
  if (right >= 0)
  {
    if (i == ix4 - 1)
    {
      f1[k * nx * ny + j * nx + i] -= C[n] * fieldCImu_right[k * ny + j] / hx / hx;
    }
  }
  if (top >= 0)
  {
    if (j == iy1)
    {
      f1[k * nx * ny + j * nx + i] -= C[n] * fieldCImu_top[k * nx + i] / hy / hy;
    }
  }
  if (bottom >= 0)
  {
    if (j == iy4 - 1)
    {
      f1[k * nx * ny + j * nx + i] -= C[n] * fieldCImu_bottom[k * nx + i] / hy / hy;
    }
  }
  if (front >= 0)
  {
    if (k == iz1)
    {
      f1[k * nx * ny + j * nx + i] -= C[n] * fieldCImu_front[j * nx + i] / hz / hz;
    }
  }
  if (back >= 0)
  {
    if (k == iz4 - 1)
    {
      f1[k * nx * ny + j * nx + i] -= C[n] * fieldCImu_back[j * nx + i] / hz / hz;
    }
  }
}

__global__ void
ac_updateU_new(int n, Dtype *fieldE, Dtype *fieldE1, Dtype *DDX, Dtype *DDY, Dtype *DDZ, Dtype LE, Dtype KE, Dtype epn2, Dtype *phiE)
{
  double tmp, Hijk;
  int i, j, k, l;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  l = k * ny * nx + j * nx + i;
  tmp = kkz * DDZ[k] + kky * DDY[j] + kkx * DDX[i];
  Hijk = -LE * (tmp * epn2 - KE);
  if (fabs(Hijk) < 1.0e-8)
  {
    Hijk = 0.0;
  }
  if (fabs(Hijk) > 1.0e-8)
  {
    tmp = exp(-dt * Hijk);
    phiE[l] = (1.0 - tmp) / Hijk;
  }
  else
  {
    phiE[l] = dt;
  }
  tmp = 1.0 - phiE[l] * Hijk;
  fieldE[l] = tmp * fieldE[l] + phiE[l] * fieldE1[l];
}

__global__ void
ch_update_phi(int n, double *phiCI, double *DDX, double *DDY, double *DDZ, double LE, double KE, double epn2, double *C)
{
  double tmp, Hijk;
  int i, j, k, l;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  double v;
  l = k * ny * nx + j * nx + i;
  tmp = DDZ[k] + DDY[j] + DDX[i];
  Hijk = LE * (tmp * tmp * epn2 - KE * tmp);
  if (n == 0)
  {
    v = C[0];
    if (fabs(Hijk) < 1.0e-8 || fabs(v) < 1.0e-8)
    {
      phiCI[l] = dt;
    }
    else
    {
      tmp = exp(-dt * Hijk * v);
      phiCI[l] = (1.0 - tmp) / (Hijk * v);
    }
  }
  else
  {
    v = C[1];
    if (fabs(Hijk) < 1.0e-8 || fabs(v) < 1.0e-8)
    {
      phiCI[l] = dt;
    }
    else
    {
      tmp = exp(-dt * Hijk * v);
      phiCI[l] = (1.0 - tmp) / (Hijk * v);
    }
  }
}

__global__ void
ch_updateU_new(int n, double *fieldCI, double *fieldCI1, double *phiCI, double *DDX, double *DDY, double *DDZ, double LE, double KE, double epn2, double *C)
{
  double tmp, Hijk;
  int i, j, k, l;
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  double v;
  l = k * ny * nx + j * nx + i;
  tmp = DDZ[k] + DDY[j] + DDX[i];
  Hijk = LE * (tmp * tmp * epn2 - KE * tmp);
  if (n == 0)
  {
    v = C[0];
    if (fabs(Hijk) < 1.0e-8 || fabs(v) < 1.0e-8)
    {
      Hijk = 0.0;
      v = 0.0;
    }
    tmp = 1.0 - phiCI[l] * Hijk * v;
    fieldCI[l] = tmp * fieldCI[l] + phiCI[l] * fieldCI1[l];
  }
  else if (n == 1)
  {
    v = C[1];
    if (fabs(Hijk) < 1.0e-8 || fabs(v) < 1.0e-8)
    {
      Hijk = 0.0;
      v = 0.0;
    }
    tmp = 1.0 - phiCI[l] * Hijk * v;
    fieldCI[l] = tmp * fieldCI[l] + phiCI[l] * fieldCI1[l];
  }
}

__global__ void
cal_new_BN(Dtype *BN1, Dtype *BN)
{
  int i, j, k;
  int num = Approx;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (myrank == 0)
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[k * NX * NY + j * NX + i];
  }
  else if (myrank == procs[0] - 1)
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[(k+(NX-Approx)) * NX * NY + j * NX + i];
  }
  else if (myrank == (procs[1] - 1) * procs[0])
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[k * NX * NY + (j+(NX-Approx)) * NX + i];
  }
  else if (myrank == (procs[1] - 1) * procs[0] + (procs[0] - 1))
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[(k+(NX-Approx)) * NX * NY + (j+(NX-Approx)) * NX + i];
  }
  else if (myrank == (procs[2] - 1) * procs[1] * procs[0])
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[k * NX * NY + j * NX + (i+(NX-Approx))];
  }
  else if (myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[0] - 1))
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[(k+(NX-Approx)) * NX * NY + j * NX + (i+(NX-Approx))];
  }
  else if (myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0])
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[k * NX * NY + (j+(NX-Approx)) * NX + (i+(NX-Approx))];
  }
  else if (myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0] + (procs[0] - 1))
  {
    BN1[k * num * num + j * num + i] = 0.0; // BN[(k+(NX-Approx)) * NX * NY + (j+(NX-Approx)) * NX + (i+(NX-Approx))];
  }
}

__global__ void
elastic_multiply_BN(Dtype *Bnre, Dtype *Bnim, Dtype *BN, Dtype *ftre, Dtype *ftim, int q)
{
  int i, j, k;
  int num = 2 * Approx;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  if (q == 0)
  {
    Bnre[k * num * num + j * num + i] = BN[k * num * num + j * num + i] * ftre[k * num * num + j * num + i];
    Bnim[k * num * num + j * num + i] = BN[k * num * num + j * num + i] * ftim[k * num * num + j * num + i];
  }
  else
  {
    Bnre[k * num * num + j * num + i] += BN[k * num * num + j * num + i] * ftre[k * num * num + j * num + i];
    Bnim[k * num * num + j * num + i] += BN[k * num * num + j * num + i] * ftim[k * num * num + j * num + i];
  }
}

__global__ void
elastic_calculate_sigma2d(Dtype *Sigma2d1, Dtype *C4D1, Dtype *Epsilon2d1)
{
  int i, j, k, l, m;
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;
  if (i < DIM && j < DIM && k == 0)
  {
    for (m = 0; m < DIM; m++)
    {
      for (l = 0; l < DIM; l++)
      {
        Sigma2d1[i * DIM + j] += C4D1[i * DIM * DIM * DIM + j * DIM * DIM + m * DIM + l] * Epsilon2d1[m * DIM + l];
      }
    }
  }
}
__global__ void
elastic_calculate_n(Dtype *n11, Dtype *n22, Dtype *n33)
{
  int i, j, k;
  int num = Approx;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  int gx, gy, gz;
  int gnnx = procs[0] * NX;
  int cntx = gnnx / 2;
  int gnny = procs[1] * NY;
  int cnty = gnny / 2;
  int gnnz = procs[2] * NZ;
  int cntz = gnnz / 2;
  if (myrank == 0 || myrank == procs[0] - 1 || myrank == (procs[1] - 1) * procs[0] || myrank == (procs[1] - 1) * procs[0] + (procs[0] - 1))
  {
    if (i < Approx && j == 0 && k == 0)
    {
      gx = cart_id[0] * NX + i;
      n11[i] = 1.0 * gx;
    }
  }
  else if (myrank == (procs[2] - 1) * procs[1] * procs[0] || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[0] - 1) || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0] || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0] + (procs[0] - 1))
  {
    if (i < Approx && j == 0 && k == 0)
    {
      gx = cart_id[0] * NX + i + (NX - num);
      n11[i] = -1.0 * (gnnx - gx);
    }
  }
  if (myrank == 0 || myrank == procs[0] - 1 || myrank == (procs[2] - 1) * procs[1] * procs[0] || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[0] - 1))
  {
    if (i == 0 && j < Approx && k == 0)
    {
      gy = cart_id[1] * NY + j;
      n22[j] = 1.0 * gy;
    }
  }
  else if (myrank == (procs[1] - 1) * procs[0] || myrank == (procs[1] - 1) * procs[0] + (procs[0] - 1) || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0] || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0] + (procs[0] - 1))
  {
    if (i == 0 && j < Approx && k == 0)
    {
      gy = cart_id[1] * NY + j + (NY - num);
      n22[j] = -1.0 * (gnny - gy);
    }
  }
  if (myrank == 0 || myrank == (procs[1] - 1) * procs[0] || myrank == (procs[2] - 1) * procs[1] * procs[0] || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0])
  {
    if (i == 0 && j == 0 && k < Approx)
    {
      gz = cart_id[2] * NZ + k;
      n33[k] = 1.0 * gz;
    }
  }
  else if (myrank == procs[0] - 1 || myrank == (procs[1] - 1) * procs[0] + (procs[0] - 1) || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[0] - 1) || myrank == (procs[2] - 1) * procs[1] * procs[0] + (procs[1] - 1) * procs[0] + (procs[0] - 1))
  {
    if (i == 0 && j == 0 && k < Approx)
    {
      gz = cart_id[2] * NZ + k + (NZ - num);
      n33[k] = -1.0 * (gnnz - gz);
    }
  }
}

__global__ void
elastic_calculate_tmp0(Dtype *tmp0, Dtype *epsilon2d, Dtype *C4D, Dtype *epsilon2d2)
{
  int i, j, k, l, m;
  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;
  if (i < DIM && j < DIM && k == 0)
  {
    tmp0[i * DIM + j] = 0.0;
    for (m = 0; m < DIM; m++)
    {
      for (l = 0; l < DIM; l++)
      {
        tmp0[i * DIM + j] += C4D[i * DIM * DIM * DIM + j * DIM * DIM + m * DIM + l] * epsilon2d[i * DIM + j] * epsilon2d2[m * DIM + l];
      }
    }
  }
  if (i < DIM && j < DIM && k == 0)
  {
    for (int stride = 1; stride < DIM * DIM; stride *= 2)
    {
      if (((i * DIM + j) % (2 * stride) == 0) && ((i * DIM + j) < (DIM * DIM - 1)))
      //                        if (((i * DIM + j) % (2 * stride) == 0))
      {
        tmp0[i * DIM + j] += tmp0[i * DIM + j + stride];
      }
    }
  }
}
#if 0
__global__ void
elastic_multiply_BN_0 (Dtype *Bnre, Dtype *Bnim, Dtype *ftre, Dtype *ftim, Dtype *tmp2, Dtype *tmp1, Dtype *tmp0, Dtype *n11, Dtype *n22, Dtype *n33, Dtype *n123, Dtype *omega_inver, Dtype *omega, Dtype *C4D, Dtype *E, Dtype *Sigma2d1, Dtype *Sigma2d2, int p, int q)
{
  int num = Approx;

  for (int z = 0; z < num; z++)
  {
	for (int y = 0; y < num; y++)
	{
		for (int x = 0; x < num; x++)
		{
  			int i, j, k;
  			i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  			j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  			k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;
  			double BN[1];
			//n123
			if (i < 1 && j == 0 && k == 0)
			{
				int gx = cart_id[0] * NX + x;
				int gy = cart_id[1] * NY + y;
				int gz = cart_id[2] * NZ + z;
				if (gx == 0 && gy == 0 && gz == 0)
				{
					n123[i+0] = 0;
					n123[i+1] = 1;
					n123[i+2] = 0;
				}
				else
				{
					double norm = 1.0 / sqrt(n11[x]*n11[x]+n22[y]*n22[y]+n33[z]*n33[z]);
					n123[i+0] = n11[x] * norm;
					n123[i+1] = n22[y] * norm;
					n123[i+2] = n33[z] * norm;
				}
			}
			//omega_inver
			if (i < DIM && j < DIM && k == 0)
			{
				if (i == j)
				{
					E[i * DIM + j] = 1.0;
				}
				else
				{
					E[i * DIM + j] = 0.0;
				}
				omega_inver[i * DIM + j] = 0.0;
				for (int m = 0; m < DIM; ++m)
				{
					for (int l = 0; l < DIM; ++l)
					{
						omega_inver[i * DIM + j] += C4D[i * DIM * DIM * DIM + m * DIM * DIM + j * DIM + l] * n123[m] * n123[l];
					}
				}
			}	
			//omega
#if 1
			if (i < DIM && j == 0 && k == 0)
			{
				double As, Bs;
				for (int n = 0; n < DIM; n++)
				{
					As = omega_inver[n * DIM + n];
					omega_inver[n * DIM + i] = omega_inver[n * DIM + i] / As;
					if (n == 0)
					{
						omega[n * DIM + i] = E[n * DIM + i] / As;
					}
					else
					{
						omega[n * DIM + i] = omega[n * DIM + i] / As;
					}
					for (int s = 0; s < DIM; s++)
					{
						if (s != n)
						{
							Bs = omega_inver[s * DIM + n];
							omega_inver[s * DIM + i] = omega_inver[s * DIM + i] - omega_inver[n * DIM + i] * Bs;
							if (n == 0)
							{
								omega[s * DIM + i] = E[s * DIM + i] - omega[n * DIM + i] * Bs;
							}
							else
							{
								omega[s * DIM + i] = omega[s * DIM + i] - omega[n * DIM + i] * Bs;
							}
						}
					}
				}
			}
#endif
			//tmp0-tmp1
			if (i < DIM && j < DIM && k == 0)
			{
				tmp1[i * DIM + j] = 0.0;
				for (int m=0; m < DIM; ++m)
				{
					for (int l = 0; l < DIM; ++l)
					{
						tmp1[i * DIM + j] += n123[i] * Sigma2d1[i * DIM + j] * omega[j * DIM + m] * Sigma2d2[m * DIM + l] * n123[l];
					}
				}
			}	
			if (i < DIM && j < DIM && k == 0)
			{
				for (int stride = 1; stride < DIM * DIM; stride *= 2)
				{
					if (((i * DIM + j) % (2 * stride) == 0) && ((i * DIM + j) < (DIM * DIM - 1)))
					//if (((i * DIM + j) % (2 * stride) == 0))
					{
						tmp1[i * DIM + j] += tmp1[i * DIM + j + stride];
					}
				}
				tmp2[i * DIM + j] = tmp0[i * DIM + j] - tmp1[i * DIM + j];
			}	
			if (myrank == 0)
			{
				if (i < 1 && j < 1 && k < 1)
				{			
					//BN[i] = tmp0[i * DIM + j] - tmp1[i * DIM + j];
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i)] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i)] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i)] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i)] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
				}
			} 
			else if (myrank == procs[0]-1)
			{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
  	  	  	else if (myrank == (procs[1]-1)*procs[0])
  	  	  	{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 2 * ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 2 * ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 2 * ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 2 * ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
  	  	  	else if (myrank == (procs[1]-1)*procs[0]+(procs[0]-1))
  	  	  	{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 3 * ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 3 * ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 3 * ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 3 * ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
  	  	  	else if (myrank == (procs[2]-1)*procs[1]*procs[0])
  	  	  	{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 4 * ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 4 * ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 4 * ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 4 * ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
  	  	  	else if (myrank == (procs[2]-1)*procs[1]*procs[0]+(procs[0]-1))
  	  	  	{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 5 * ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 5 * ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 5 * ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 5 * ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
  	  	  	else if (myrank == (procs[2]-1)*procs[1]*procs[0]+(procs[1]-1)*procs[0])
  	  	  	{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 6 * ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 6 * ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 6 * ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 6 * ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
  	  	  	else if (myrank == (procs[2]-1)*procs[1]*procs[0]+(procs[1]-1)*procs[0]+(procs[0]-1))
  	  	  	{
				if (i < 1 && j < 1 && k < 1)
				{			
//					B[(k+z)*num*num+(j+y)*num+x+i] = tmp2[0];
					if (q == 0) {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 7 * ela_size] = tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 7 * ela_size] = tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					}
					else {
						Bnre[(z+k) * num * num + (y+j) * num + (x+i) + 7 * ela_size] += tmp2[0] * ftre[(z+k) * num * num + (y+j) * num + (x+i)];
						Bnim[(z+k) * num * num + (y+j) * num + (x+i) + 7 * ela_size] += tmp2[0] * ftim[(z+k) * num * num + (y+j) * num + (x+i)];
					} 
				}
			}
		}
	}
  }
}
#endif

__global__ void
elastic_copyin(Dtype *f, Dtype *fe)
{
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f[k * NX * NY + j * NX + i] = fe[(k + nghost) * ny * nx + (j + nghost) * nx + i + nghost];
}

__global__ void
elastic_copyout(Dtype *f, Dtype *fe)
{
  int i, j, k;
  double np0 = (double)(NX * procs[0]);
  double np1 = (double)(NY * procs[1]);
  double np2 = (double)(NZ * procs[2]);
  double c = 1.0 / np0 / np1 / np2;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f[(k + nghost) * ny * nx + (j + nghost) * nx + i + nghost] += c * fe[k * NX * NY + j * NX + i];
  // f[(k + nghost) * ny * nx + (j + nghost) * nx + i + nghost] += fe[k * NX * NY + j * NX + i];
}

__global__ void
down_sample_eta(Dtype *f, Dtype *fe)
{
  int i, j, k;
  int dim = (nx - 2 * nghost) / 1;
  int stride = 1;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f[k * dim * dim + j * dim + i] = fe[(k * stride + nghost) * ny * nx + (j * stride + nghost) * nx + i * stride + nghost];
}

__global__ void
initial_fftre_im(double *re_out, double *im_out)
{
  int i, j, k;
  int num = Approx * 2;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  re_out[k * num * num + j * num + i] = 0.0;
  im_out[k * num * num + j * num + i] = 0.0;
}

__global__ void
IGx_IGy_IGz_matrix(double *Gx_re, double *Gx_im, double *Gy_re, double *Gy_im, double *Gz_re, double *Gz_im, int flag)
{
  int i, j, k;
  int dim = (nx - 2 * nghost) / 1;
  int num = Approx;

  double dc0 = (double)dim * cart_id[0];
  double dc1 = (double)dim * cart_id[1];
  double dc2 = (double)dim * cart_id[2];
  double dp0 = (double)dim * procs[0];
  double dp1 = (double)dim * procs[1];
  double dp2 = (double)dim * procs[2];

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;

  if (flag == 0)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 4)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 1)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
  else if (flag == 5)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
  else if (flag == 2)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 6)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 3)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
  else if (flag == 7)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
}

__global__ void
Gx_Gy_Gz_matrix(double *Gx_re, double *Gx_im, double *Gy_re, double *Gy_im, double *Gz_re, double *Gz_im, int flag)
{
  int i, j, k;
  int dim = (nx - 2 * nghost) / 1;
  int num = Approx;

  double dc0 = (double)dim * cart_id[0];
  double dc1 = (double)dim * cart_id[1];
  double dc2 = (double)dim * cart_id[2];
  double dp0 = (double)dim * procs[0];
  double dp1 = (double)dim * procs[1];
  double dp2 = (double)dim * procs[2];

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;

  if (flag == 0)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 4)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 1)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
  else if (flag == 5)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * j / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * j / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
  else if (flag == 2)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 6)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * j / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * j / dp2);
  }
  else if (flag == 3)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * j / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * j / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
  else if (flag == 7)
  {
    Gx_re[j * dim + i] = cos(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gx_im[j * dim + i] = sin(2.0 * PI * (dc0 + i) * (dp0 + j - num) / dp0);
    Gy_re[j * dim + i] = cos(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gy_im[j * dim + i] = sin(2.0 * PI * (dc1 + i) * (dp1 + j - num) / dp1);
    Gz_re[j * dim + i] = cos(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
    Gz_im[j * dim + i] = sin(2.0 * PI * (dc2 + i) * (dp2 + j - num) / dp2);
  }
}

__global__ void
initial_elas(double *elas)
{
  int i, j, k;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  elas[k * nx * ny + j * nx + i] = 0.0;
}

void elastic_calculate_forward(int p)
{
  dim3 blocks((ix4 - ix1) / THREADS_PER_BLOCK_X, (iy4 - iy1) / THREADS_PER_BLOCK_Y, (iz4 - iz1) / THREADS_PER_BLOCK_Z);
  dim3 threads(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela1((NX / 1) / THREADS_PER_BLOCK_X, (NY / 1) / THREADS_PER_BLOCK_Y, (NZ / 1) / THREADS_PER_BLOCK_Z);
  dim3 threads_ela1(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela2(NX / THREADS_PER_BLOCK_X, NY / THREADS_PER_BLOCK_Y, NZ / THREADS_PER_BLOCK_Z);
  dim3 threads_ela2(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela(NX / THREADS_PER_BLOCK_X, Approx / THREADS_PER_BLOCK_Y, Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_ela(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela_16(Approx / THREADS_PER_BLOCK_X, Approx / THREADS_PER_BLOCK_Y, Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_ela_16(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela3(2 * Approx / THREADS_PER_BLOCK_X, 2 * Approx / THREADS_PER_BLOCK_Y, 2 * Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_ela3(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);

  int flag;
  hipLaunchKernelGGL(initial_fftre_im, blocks_ela3, threads_ela3, 0, 0, tmpy_fftRE[p], tmpy_fftIM[p]);
  hipLaunchKernelGGL(down_sample_eta, blocks_ela1, threads_ela1, 0, 0, tmpy_RE1, fieldE + p * offset);
  for (flag = 0; flag < 8; flag++)
  {
    hipLaunchKernelGGL(IGx_IGy_IGz_matrix, blocks_ela, threads_ela, 0, 0, Gx_re, Gx_im, Gy_re, Gy_im, Gz_re, Gz_im, flag);
    // x -- dim
    //(tmpy_RE1)^T*(Gx_re)^T <--> ((ny*nz)*nx)*(nx*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gx_n, Gx_k, Gx_m, &alpha, tmpy_RE1, Gx_m, Gx_re, Gx_m, &beta, temp1, Gx_n);
    //(tmpy_RE1)^T*(Gx_im)^T <--> ((ny*nz)*nx)*(nx*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gx_n, Gx_k, Gx_m, &alpha1, tmpy_RE1, Gx_m, Gx_im, Gx_m, &beta, temp2, Gx_n);
    // y -- dim
    //(temp1)^T*(Gy_re)^T <--> ((16*nz)*ny)*(ny*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gy_n, Gy_k, Gy_m, &alpha, temp1, Gy_m, Gy_re, Gy_m, &beta, temp3, Gy_n);
    //(temp2)^T*(Gy_im)^T <--> ((16*nz)*ny)*(ny*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gy_n, Gy_k, Gy_m, &alpha, temp2, Gy_m, Gy_im, Gy_m, &beta, temp4, Gy_n);
    // caculate real part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gy_n, Gy_k, &alpha, temp3, Gy_n, &beta2, temp4, Gy_n, temp5, Gy_n);

    //(temp1)^T*(Gy_im)^T <--> ((16*nz)*ny)*(ny*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gy_n, Gy_k, Gy_m, &alpha1, temp1, Gy_m, Gy_im, Gy_m, &beta, temp3, Gy_n);
    //(temp2)^T*(Gy_re)^T <--> ((16*nz)*ny)*(ny*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gy_n, Gy_k, Gy_m, &alpha, temp2, Gy_m, Gy_re, Gy_m, &beta, temp4, Gy_n);
    // caculate imag part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gy_n, Gy_k, &alpha, temp3, Gy_n, &beta2, temp4, Gy_n, temp6, Gy_n);

    // z -- dim
    //(temp7)^T*(Gz_re)^T <--> ((16*16)*nz)*(nz*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gz_n, Gz_k, Gz_m, &alpha, temp5, Gz_m, Gz_re, Gz_m, &beta, temp7, Gz_n);
    //(temp8)^T*(Gz_im)^T <--> ((16*16)*nz)*(nz*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gz_n, Gz_k, Gz_m, &alpha, temp6, Gz_m, Gz_im, Gz_m, &beta, temp8, Gz_n);
    // caculate real part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gz_n, Gz_k, &alpha, temp7, Gz_n, &beta2, temp8, Gz_n, tmpy_fftRE[p] + flag * ela_size, Gz_n);

    //(temp7)^T*(Gz_im)^T <--> ((16*16)*nz)*(nz*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gz_n, Gz_k, Gz_m, &alpha1, temp5, Gz_m, Gz_im, Gz_m, &beta, temp7, Gz_n);
    //(temp8)^T*(Gz_re)^T <--> ((16*16)*nz)*(nz*16)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, Gz_n, Gz_k, Gz_m, &alpha, temp6, Gz_m, Gz_re, Gz_m, &beta, temp8, Gz_n);
    // caculate imag part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gz_n, Gz_k, &alpha, temp7, Gz_n, &beta2, temp8, Gz_n, tmpy_fftIM[p] + flag * ela_size, Gz_n);
  }

  hipMemcpy(tmpy_fftre[p], tmpy_fftRE[p], sizeof(Dtype) * 8 * ela_size, hipMemcpyDeviceToHost);
  hipMemcpy(tmpy_fftim[p], tmpy_fftIM[p], sizeof(Dtype) * 8 * ela_size, hipMemcpyDeviceToHost);

  MPI_Allreduce(&tmpy_fftre[p][0], &fftre[p][0], 8 * ela_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&tmpy_fftim[p][0], &fftim[p][0], 8 * ela_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  hipMemcpy(fftRE[p], fftre[p], sizeof(Dtype) * 8 * ela_size, hipMemcpyHostToDevice);
  hipMemcpy(fftIM[p], fftim[p], sizeof(Dtype) * 8 * ela_size, hipMemcpyHostToDevice);
}
void elastic_calculate_backward(int p)
{
  dim3 blocks((ix4 - ix1) / THREADS_PER_BLOCK_X, (iy4 - iy1) / THREADS_PER_BLOCK_Y, (iz4 - iz1) / THREADS_PER_BLOCK_Z);
  dim3 threads(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela1((NX / 1) / THREADS_PER_BLOCK_X, (NY / 1) / THREADS_PER_BLOCK_Y, (NZ / 1) / THREADS_PER_BLOCK_Z);
  dim3 threads_ela1(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela2(NX / THREADS_PER_BLOCK_X, NY / THREADS_PER_BLOCK_Y, NZ / THREADS_PER_BLOCK_Z);
  dim3 threads_ela2(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela(NX / THREADS_PER_BLOCK_X, Approx / THREADS_PER_BLOCK_Y, Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_ela(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela_16(Approx / THREADS_PER_BLOCK_X, Approx / THREADS_PER_BLOCK_Y, Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_ela_16(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks_ela3(2 * Approx / THREADS_PER_BLOCK_X, 2 * Approx / THREADS_PER_BLOCK_Y, 2 * Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_ela3(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  for (int q = 0; q < NELASTIC; q++)
  {
    hipLaunchKernelGGL(elastic_multiply_BN, blocks_ela3, threads_ela3, 0, 0, BN_fftRE[p], BN_fftIM[p], BN1[p][q], fftRE[q], fftIM[q], q);
  }
  hipLaunchKernelGGL(initial_elas, blocks, threads, 0, 0, Elas);
  for (flag = 0; flag < 8; flag++)
  {
    hipLaunchKernelGGL(Gx_Gy_Gz_matrix, blocks_ela, threads_ela, 0, 0, Gx_re, Gx_im, Gy_re, Gy_im, Gz_re, Gz_im, flag);
    // x -- dim
    //(BN_fftRE[p])^T*(Gx_re) <--> ((16*16)*16)*(16*nx)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gz_n, Gz_m, Gz_k, &alpha, fftRE[p] + flag * ela_size, Gz_k, Gx_re, Gz_m, &beta, temp3, Gz_n);
    //(BN_fftIM[p])^T*(Gx_im) <--> ((16*16)*16)*(16*nx)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gz_n, Gz_m, Gz_k, &alpha1, fftIM[p] + flag * ela_size, Gz_k, Gx_im, Gz_m, &beta, temp4, Gz_n);
    // caculate real part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gz_n, Gz_m, &alpha, temp3, Gz_n, &beta2, temp4, Gz_n, temp5, Gz_n);

    //(BN_fftRE[p])^T*(Gx_im) <--> ((16*16)*16)*(16*nx)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gz_n, Gz_m, Gz_k, &alpha, fftRE[p] + flag * ela_size, Gz_k, Gx_im, Gz_m, &beta, temp3, Gz_n);
    //(BN_fftIM[p])^T*(Gx_re) <--> ((16*16)*16)*(16*nx)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gz_n, Gz_m, Gz_k, &alpha, fftIM[p] + flag * ela_size, Gz_k, Gx_re, Gz_m, &beta, temp4, Gz_n);
    // caculate imag part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gz_n, Gz_m, &alpha, temp3, Gz_n, &beta2, temp4, Gz_n, temp6, Gz_n);

    // y -- dim
    //(temp3)^T*(Gy_re) <--> ((16*nx)*16)*(16*ny)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gy_n, Gy_m, Gy_k, &alpha, temp5, Gy_k, Gy_re, Gy_m, &beta, temp1, Gy_n);
    //(temp4)^T*(Gy_im) <--> ((16*nx)*16)*(16*ny)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gy_n, Gy_m, Gy_k, &alpha1, temp6, Gy_k, Gy_im, Gy_m, &beta, temp2, Gy_n);
    // caculate real part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gy_n, Gy_m, &alpha, temp1, Gy_n, &beta2, temp2, Gy_n, Itemp1, Gy_n);
    //(temp3)^T*(Gy_im) <--> ((16*nx)*16)*(16*ny)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gy_n, Gy_m, Gy_k, &alpha, temp5, Gy_k, Gy_im, Gy_m, &beta, temp1, Gy_n);
    //(temp4)^T*(Gy_re) <--> ((16*nx)*16)*(16*ny)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gy_n, Gy_m, Gy_k, &alpha, temp6, Gy_k, Gy_re, Gy_m, &beta, temp2, Gy_n);
    // caculate imag part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gy_n, Gy_m, &alpha, temp1, Gy_n, &beta2, temp2, Gy_n, Itemp2, Gy_n);
    // z -- dim
    //(Itemp1)^T*(Gz_re) <--> ((nx*ny)*16)*(16*nz)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gx_n, Gx_m, Gx_k, &alpha, Itemp1, Gx_k, Gz_re, Gx_m, &beta, Itemp3, Gx_n);
    //(Itemp2)^T*(Gz_im) <--> ((nx*ny)*16)*(16*nz)
    rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, Gx_n, Gx_m, Gx_k, &alpha1, Itemp2, Gx_k, Gz_im, Gx_m, &beta, Itemp4, Gx_n);
    // caculate real part
    rocblas_dgeam(handle, rocblas_operation_none, rocblas_operation_none, Gx_n, Gx_m, &alpha, Itemp3, Gx_n, &beta2, Itemp4, Gx_n, tmpy_RE1, Gx_n);
    hipLaunchKernelGGL(elastic_copyout, blocks_ela2, threads_ela2, 0, 0, Elas, tmpy_RE1);
  }
}

void elastic_transfer(int n)
{
  int threads_x, threads_y, threads_z;
  //  for (int n = 0; n < nac; n++)
  //  {
  hipMemcpy(s_front, &Elas[nx * ny * nghost], sizeof(Dtype) * fb_size, hipMemcpyDeviceToHost);
  hipMemcpy(s_back, &Elas[nx * ny * (nz - nghost - (nghost + 2))], sizeof(Dtype) * fb_size, hipMemcpyDeviceToHost);

  MPI_Startall(4, ireq_front_back);
  MPI_Waitall(4, ireq_front_back, status);

  hipMemcpy(R_front, r_front, sizeof(Dtype) * fb_size, hipMemcpyHostToDevice);
  hipMemcpy(R_back, r_back, sizeof(Dtype) * fb_size, hipMemcpyHostToDevice);
  dim3 blocks_tb_pack(nx / THREADS_PER_BLOCK_X, (nghost + 2) / THREADS_PER_BLOCK_Y, (nz + (nghost + 2) * 2) / THREADS_PER_BLOCK_Z);
  dim3 threads_tb_pack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  hipLaunchKernelGGL(top_bottom_pack, blocks_tb_pack, threads_tb_pack, 0, 0,
                     Elas, S_top, S_bottom, R_front, R_back);

  hipMemcpy(s_top, S_top, sizeof(Dtype) * tb_size, hipMemcpyDeviceToHost);
  hipMemcpy(s_bottom, S_bottom, sizeof(Dtype) * tb_size, hipMemcpyDeviceToHost);
  MPI_Startall(4, ireq_top_bottom);
  MPI_Waitall(4, ireq_top_bottom, status);

  hipMemcpy(R_top, r_top, sizeof(Dtype) * tb_size, hipMemcpyHostToDevice);
  hipMemcpy(R_bottom, r_bottom, sizeof(Dtype) * tb_size, hipMemcpyHostToDevice);

  dim3 blocks_lr_pack((nghost + 2) / THREADS_PER_BLOCK_X, (ny + (nghost + 2) * 2) / THREADS_PER_BLOCK_Y, (nz + (nghost + 2) * 2) / THREADS_PER_BLOCK_Z);
  dim3 threads_lr_pack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  hipLaunchKernelGGL(left_right_pack, blocks_lr_pack, threads_lr_pack, 0, 0,
                     Elas, S_left, S_right, R_top, R_bottom, R_front, R_back);

  hipMemcpy(s_left, S_left, sizeof(Dtype) * lr_size, hipMemcpyDeviceToHost);
  hipMemcpy(s_right, S_right, sizeof(Dtype) * lr_size, hipMemcpyDeviceToHost);

  MPI_Startall(4, ireq_left_right);
  MPI_Waitall(4, ireq_left_right, status);

  hipMemcpy(R_left, r_left, sizeof(Dtype) * lr_size, hipMemcpyHostToDevice);
  hipMemcpy(R_right, r_right, sizeof(Dtype) * lr_size, hipMemcpyHostToDevice);

  threads_x = nghost;
  dim3 blocks_lr_unpack(nghost / threads_x, ny / THREADS_PER_BLOCK_Y, nz / THREADS_PER_BLOCK_Z);
  dim3 threads_lr_unpack(threads_x, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  hipLaunchKernelGGL(left_right_unpack, blocks_lr_unpack, threads_lr_unpack, 0, 0,
                     Elas, R_left, R_right);
  threads_y = nghost;
  dim3 blocks_tb_unpack(nx / THREADS_PER_BLOCK_X, nghost / threads_y, nz / THREADS_PER_BLOCK_Z);
  dim3 threads_tb_unpack(THREADS_PER_BLOCK_X, threads_y, THREADS_PER_BLOCK_Z);
  hipLaunchKernelGGL(top_bottom_unpack, blocks_tb_unpack, threads_tb_unpack, 0, 0,
                     Elas, R_top, R_bottom);
  threads_z = nghost;
  dim3 blocks_fb_unpack(nx / THREADS_PER_BLOCK_X, ny / THREADS_PER_BLOCK_Y, nghost / threads_z);
  dim3 threads_fb_unpack(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, threads_z);
  hipLaunchKernelGGL(front_back_unpack, blocks_fb_unpack, threads_fb_unpack, 0, 0,
                     Elas, R_front, R_back);
  //  }
}

__global__ void
skip_copyin_new_2(Dtype *f, Dtype *fe)
{
  int i, j, k;
  int dim_x = (nx - nghost * 2) / 2;
  int dim_y = (ny - nghost * 2) / 2;
  int stride = 2;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  f[k * dim_y * dim_x + j * dim_x + i] = (1.0 / 8.0) *
                                         (fe[((k * stride) + nghost) * ny * nx + ((j * stride) + nghost) * nx + (i * stride) + nghost] +
                                          fe[((k * stride) + nghost) * ny * nx + ((j * stride) + nghost) * nx + (i * stride + 1) + nghost] +
                                          fe[((k * stride) + nghost) * ny * nx + ((j * stride + 1) + nghost) * nx + (i * stride) + nghost] +
                                          fe[((k * stride) + nghost) * ny * nx + ((j * stride + 1) + nghost) * nx + (i * stride + 1) + nghost] +
                                          fe[((k * stride + 1) + nghost) * ny * nx + ((j * stride) + nghost) * nx + (i * stride) + nghost] +
                                          fe[((k * stride + 1) + nghost) * ny * nx + ((j * stride) + nghost) * nx + (i * stride + 1) + nghost] +
                                          fe[((k * stride + 1) + nghost) * ny * nx + ((j * stride + 1) + nghost) * nx + (i * stride) + nghost] +
                                          fe[((k * stride + 1) + nghost) * ny * nx + ((j * stride + 1) + nghost) * nx + (i * stride + 1) + nghost]);
}

__global__ void
Dtype2B_E(Dtype *E, Stype *F2B)
{
  int i, j, k, l;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  l = k * nx * ny + j * nx + i;
  F2B[l] = (Stype)((E[l] + 0.1) * 30000);
}

__global__ void
Dtype2B_CI(Dtype *CI, Stype *F2B)
{
  int i, j, k, l;

  i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  j = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
  k = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;

  l = k * nx * ny + j * nx + i;
  F2B[l] = (Stype)(CI[l] * 50000);
}
extern void hip_init()
{
  deviceId = myrank % 4;
  // hipGetDevice(&deviceId);
  hipGetDeviceProperties(&props, deviceId);
}

// f1
void ac_calc_F1(float *field2_all, int out_iter, int detect_iter)
{

  int threads_x, threads_y, threads_z;
  dim3 blocks((ix4 - ix1) / THREADS_PER_BLOCK_X, (iy4 - iy1) / THREADS_PER_BLOCK_Y, (iz4 - iz1) / THREADS_PER_BLOCK_Z);
  dim3 threads(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);

  dim3 blocks2B(nx / THREADS_PER_BLOCK_X, ny / THREADS_PER_BLOCK_Y, nz / THREADS_PER_BLOCK_Z);
  dim3 threads2B(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);

  dim3 blocks_16(Approx / THREADS_PER_BLOCK_X, Approx / THREADS_PER_BLOCK_Y, Approx / THREADS_PER_BLOCK_Z);
  dim3 threads_16(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);
  dim3 blocks2B2(((nx - nghost * 2) / 2) / THREADS_PER_BLOCK_X, ((ny - nghost * 2) / 2) / THREADS_PER_BLOCK_Y, ((nz - nghost * 2) / 2) / THREADS_PER_BLOCK_Z);
  dim3 threads2B2(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, THREADS_PER_BLOCK_Z);

  // double tt;
  double functime, runtime, walltime;
  //int ot2 = 0;
  //int ot1 = 0;
  //iter = 0;
  //tt = 0.0;
  //	flag = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  walltime = MPI_Wtime();
  int n;
  for (n = 0; n < nac; n++)
  {
    hipMemcpy(lambda + n * 3 * 3, ac[n].lambda, 3 * 3 * sizeof(Dtype), hipMemcpyHostToDevice);
  }

  // CPU --> DCU dgemm
  hipMemcpy(mpxi, MPXI, nx * nx * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(mpyi, MPYI, ny * ny * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(mpzi, MPZI, nz * nz * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(mpx, MPX, nx * nx * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(mpy, MPY, ny * ny * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(mpz, MPZ, nz * nz * sizeof(Dtype), hipMemcpyHostToDevice);
  // CPU --> DCU update
  hipMemcpy(ddx, DDX, nx * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(ddy, DDY, ny * sizeof(Dtype), hipMemcpyHostToDevice);
  hipMemcpy(ddz, DDZ, nz * sizeof(Dtype), hipMemcpyHostToDevice);

  //while (tt < 10)
  int gyq = 0;
  while (gyq < 100)
  { // here
    dgemm_time = 0.0f;
    F1_time = 0.0f;
    F2_time = 0.0f;
    trans_time = 0.0f;
    trans_Memcpy_time = 0.0f;
    trans_MPI_time = 0.0f;
    trans_pack_time = 0.0f;
    trans_enlarge_time = 0.0f;
    trans_mu_time = 0.0f;
    trans_unpack_time = 0.0f;
    updateU_new_time = 0.0f;


    if (iter > 0)
    {
      for (int m = 0; m < nac; m++)
      {
        hipMemcpy(ac[m].fieldE, fieldE + m * offset, nx * ny * nz * sizeof(double), hipMemcpyDeviceToHost);
      }
      #if 0
      for (int m = 0; m < nch; m++)
      {
        hipMemcpy(ch[m].fieldCI, fieldCI + m * offset, nx * ny * nz * sizeof(double), hipMemcpyDeviceToHost);
      }
      #endif
      check_soln_new(tt);
    }
    if (checkpoint_ == 1 && iter < 1000)
    { count_gyq++; 
      initial_ac (ac[P].fieldE, ac[Q].fieldE);
      //for (int m = 0; m < nac; m++)
      //{
        if (P == Q)
        {
            hipMemcpy(fieldE + P * offset, ac[P].fieldE, nx * ny * nz * sizeof(double), hipMemcpyHostToDevice);
        }
        else
        {
            hipMemcpy(fieldE + P * offset, ac[P].fieldE, nx * ny * nz * sizeof(double), hipMemcpyHostToDevice);
            hipMemcpy(fieldE + Q * offset, ac[Q].fieldE, nx * ny * nz * sizeof(double), hipMemcpyHostToDevice);
        }
      //}
      
    } 
    if (iter == 0)
    {
      check_soln_new(tt);
    }

    // adaptive time step

    #if 0
    for (int m = 0; m < nac; m++)
    {
      hipMemcpy(ac[m].fieldE_old, fieldE + m * offset, nx * ny * nz * sizeof(double), hipMemcpyDeviceToHost);
    }
    for (int m = 0; m < nch; m++)
    {
      hipMemcpy(ch[m].fieldCI_old, fieldCI + m * offset, nx * ny * nz * sizeof(double), hipMemcpyDeviceToHost);
    }
    #endif

    // adaptive time step

    tt += dt;
    runtime = MPI_Wtime();
#ifdef SCLETD_DEBUG
    hipEventRecord(st2, NULL);
    hipEventSynchronize(st2);
#endif

    transfer();

#ifdef SCLETD_DEBUG
    hipEventRecord(ed2, NULL);
    hipEventSynchronize(ed2);
    hipEventElapsedTime(&timer2, st2, ed2);
    trans_time = timer2;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    for (n = 0; n < nac; n++)
    {
      elastic_calculate_forward(n);
    }
    for (n = 0; n < nac; n++)
    {
      elastic_calculate_backward(n);
      elastic_transfer(n);

      hipLaunchKernelGGL(anisotropic_calc_dev_2, blocks, threads, 0, 0,
                         f_aniso, fieldEr + n * offset_Er, lambda, Ori, n);

      hipLaunchKernelGGL(ac_calc_F1_dev, blocks, threads, 0, 0,
                         fieldE1 + n * offset, fieldE, fieldCI, fieldEu_left + n * u_lr, fieldEu_right + n * u_lr,\ 
		fieldEu_top + n * u_tb,
                         fieldEu_bottom + n * u_tb, fieldEu_front + n * u_fb,
                         fieldEu_back + n * u_fb, ac[n].LE, ac[n].KE, ac[n].epn2, n);
     //printf("Gbcc = %lf, Ghcp = %lf\n", Gbcc, Ghcp);

      hipLaunchKernelGGL(ac_calc_F2_dev, blocks, threads, 0, 0,
                         fieldE1 + n * offset, Elas, f_aniso, ac[n].LE, ac[n].epn2, F, fieldE + n * offset, n);
    }
#if 1

#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    F1_time += timer;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif
#if 0
    for (n = 0; n < nch; n++)
    {
      hipLaunchKernelGGL(update_M_C, blocks, threads, 0, 0,
                         n, C, M, fieldE, fieldCI);

      hipLaunchKernelGGL(ch_get_df, blocks, threads, 0, 0,
                         n, dfc1, dfc2, fieldE, fieldCI);

      hipLaunchKernelGGL(ch_calc_F1_A, blocks, threads, 0, 0,
                         n, ft, dfc1, dfc2, fieldCI, ch[n].KE);

      hipLaunchKernelGGL(ch_calc_F1_B, blocks, threads, 0, 0,
                         n, f1, ft, M, C);

      if (n == 0)
      {
        hipLaunchKernelGGL(ch_calc_F2, blocks, threads, 0, 0,
                           n, f2, M, dfc1);
        hipLaunchKernelGGL(ch_calc_F2, blocks, threads, 0, 0,
                           n, f3, M + 1 * offset, dfc2);
      }
      else
      {
        hipLaunchKernelGGL(ch_calc_F2, blocks, threads, 0, 0,
                           n, f2, M + 2 * offset, dfc1);
        hipLaunchKernelGGL(ch_calc_F2, blocks, threads, 0, 0,
                           n, f3, M + 3 * offset, dfc2);
      }

      hipLaunchKernelGGL(ac_add_F1_F2_dev, blocks, threads, 0, 0,
                         n, fieldCI1 + n * offset, f1, f2, f3, fieldCImu_left + n * u_lr, fieldCImu_right + n * u_lr,
                         fieldCImu_top + n * u_tb, fieldCImu_bottom + n * u_tb,
                         fieldCImu_front + n * u_fb, fieldCImu_back + n * u_fb, C);
    }
#endif
#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    F2_time += timer;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    for (n = 0; n < nac; n++)
    {
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, x_n, x_k, x_m, &alpha, fieldE + n * offset, x_k, mpxi, x_m, &beta, fieldEt, x_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, y_n, y_k, y_m, &alpha, fieldEt, y_k, mpyi, y_m, &beta, fieldEp, y_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, z_n, z_k, z_m, &alpha, fieldEp, z_k, mpzi, z_m, &beta, fieldE + n * offset, z_n);

      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, x_n, x_k, x_m, &alpha, fieldE1 + n * offset, x_k, mpxi, x_m, &beta, fieldEt, x_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, y_n, y_k, y_m, &alpha, fieldEt, y_k, mpyi, y_m, &beta, fieldEp, y_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, z_n, z_k, z_m, &alpha, fieldEp, z_k, mpzi, z_m, &beta, fieldE1 + n * offset, z_n);
    }
    #if 0
    for (n = 0; n < nch; n++)
    {
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, x_n, x_k, x_m, &alpha, fieldCI + n * offset, x_k, mpxi, x_m, &beta, fieldCIt, x_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, y_n, y_k, y_m, &alpha, fieldCIt, y_k, mpyi, y_m, &beta, fieldCIp, y_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, z_n, z_k, z_m, &alpha, fieldCIp, z_k, mpzi, z_m, &beta, fieldCI + n * offset, z_n);

      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, x_n, x_k, x_m, &alpha, fieldCI1 + n * offset, x_k, mpxi, x_m, &beta, fieldCIt, x_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, y_n, y_k, y_m, &alpha, fieldCIt, y_k, mpyi, y_m, &beta, fieldCIp, y_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, z_n, z_k, z_m, &alpha, fieldCIp, z_k, mpzi, z_m, &beta, fieldCI1 + n * offset, z_n);
    }
    #endif

#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    dgemm_time += timer;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    for (n = 0; n < nac; n++)
    {
      hipLaunchKernelGGL(ac_updateU_new, blocks, threads, 0, 0,
                         n, fieldE + n * offset, fieldE1 + n * offset, ddx, ddy, ddz, ac[n].LE, ac[n].KE, ac[n].epn2, phiE);
    }
    #if 0
    for (n = 0; n < nch; n++)
    {
      hipLaunchKernelGGL(ch_update_phi, blocks, threads, 0, 0,
                         n, phiCI, ddx, ddy, ddz, ch[n].LE, ch[n].KE, ch[n].epn2, C);

      hipLaunchKernelGGL(ch_updateU_new, blocks, threads, 0, 0,
                         n, fieldCI + n * offset, fieldCI1 + n * offset, phiCI, ddx,
                         ddy, ddz, ch[n].LE, ch[n].KE, ch[n].epn2, C);
    }
    #endif

#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    updateU_new_time += timer;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    for (n = 0; n < nac; n++)
    {
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, x_n, x_k, x_m, &alpha, fieldE + n * offset, x_k, mpx, x_m, &beta, fieldEt, x_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, y_n, y_k, y_m, &alpha, fieldEt, y_k, mpy, y_m, &beta, fieldEp, y_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, z_n, z_k, z_m, &alpha, fieldEp, z_k, mpz, z_m, &beta, fieldE + n * offset, z_n);
    }
    #if 0
    for (n = 0; n < nch; n++)
    {
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, x_n, x_k, x_m, &alpha, fieldCI + n * offset, x_k, mpx, x_m, &beta, fieldCIt, x_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, y_n, y_k, y_m, &alpha, fieldCIt, y_k, mpy, y_m, &beta, fieldCIp, y_n);
      rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_transpose, z_n, z_k, z_m, &alpha, fieldCIp, z_k, mpz, z_m, &beta, fieldCI + n * offset, z_n);
    }
    #endif
    // adaptive time step
    #if 0
    for (int m = 0; m < nac; m++)
    {
      hipMemcpy(ac[m].fieldE, fieldE + m * offset, nx * ny * nz * sizeof(double), hipMemcpyDeviceToHost);
    }
    for (int m = 0; m < nch; m++)
    {
      hipMemcpy(ch[m].fieldCI, fieldCI + m * offset, nx * ny * nz * sizeof(double), hipMemcpyDeviceToHost);
    }
    flag1 = 0;
    check_cahn_hilliard_min(flag1);
    if (flag1 == 1)
    {
      dt = dt * 0.1;
      for (int m = 0; m < nac; m++)
      {
        hipMemcpy(fieldE + m * offset, ac[m].fieldE_old, nx * ny * nz * sizeof(double), hipMemcpyHostToDevice);
      }
      for (int m = 0; m < nch; m++)
      {
        hipMemcpy(fieldCI + m * offset, ch[m].fieldCI_old, nx * ny * nz * sizeof(double), hipMemcpyHostToDevice);
      }
    }
    else if (flag1 == 2)
    {
      dt = dt * 1.1;
    }
    if (dt > 1.2)
      dt = 1.0;
    if (dt < 0.00005)
      dt = 0.00005;
    #endif
    // adaptive time step
    if (out_iter % detect_iter == 0)
    {
    for (int m = 0; m < nac; m++)
    {
      // skip_copyin_new_2 (field1, fieldE + m * offset);
      hipLaunchKernelGGL(skip_copyin_new_2, blocks2B2, threads2B2, 0, 0,
                         field1, fieldE + m * offset);
      hipMemcpy(field2, field1, ((nx - 2 * nghost) * (ny - 2 * nghost) * (nz - 2 * nghost) / 8) * sizeof(double), hipMemcpyDeviceToHost);
      for (int temp_i = 0; temp_i < 128 * 128 * 128; temp_i++)
      {
        (field2[temp_i] > 0.1) ? field2_all[temp_i + m * 128 * 128 * 128] = 1.0 : field2_all[temp_i + m * 128 * 128 * 128] = 0.0;
        // field2_all[temp_i + m * 128 * 128 * 128] = (float)field2[temp_i];
      }
    }
    }

#ifdef SCLETD_DEBUG
    hipEventRecord(ed, NULL);
    hipEventSynchronize(ed);
    hipEventElapsedTime(&timer, st, ed);
    dgemm_time += timer;
    hipEventRecord(st, NULL);
    hipEventSynchronize(st);
#endif

    ot1 += 1;
    if (ot1 >= nout)
    {
      ot1 = 0;
      irun += 1;
      for (n = 0; n < nac; n++)
      {
        	hipLaunchKernelGGL (Dtype2B_E, blocks2B, threads2B, 0, 0,\
			fieldE + n * offset, field2BE + n * offset);
      }
       write_field2B (irun, 0);
       #if 0
      for (n = 0; n < nch; n++)
      {
        	hipLaunchKernelGGL (Dtype2B_CI, blocks2B, threads2B, 0, 0,\
			fieldCI + n * offset, field2BCI + n * offset);
      }
       write_field2B (irun, 1);
       #endif
    }

    ot2 += 1;
    if (ot2 >= nchk)
    {
      ot2 = 0;
      write_chk();
      chk = (chk + 1) % 2;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (myrank == prank)
    {
#ifdef SCLETD_DEBUG
      printf("*******************************\n");
      printf("dgemm_time\t%lf\n", dgemm_time);
      printf("Allen_Cahn_nonlinear_time\t%lf\n", F1_time);
      printf("Cahn_Hilliard_nonlinear_time\t%lf\n", F2_time);
      printf("trans_time\t%lf\n", trans_time);
      printf("trans_Memcpy_time\t%lf\n", trans_Memcpy_time);
      printf("trans_MPI_time\t%lf\n", trans_MPI_time);
      printf("trans_pack_time\t%lf\n", trans_pack_time);
      printf("trans_enlarge_time\t%lf\n", trans_enlarge_time);
      printf("trans_mu_time\t%lf\n", trans_mu_time);
      printf("trans_unpack_time\t%lf\n", trans_unpack_time);
      printf("updateU_new_time\t%lf\n", updateU_new_time);
#endif
      printf("runtime\t%lf\n", MPI_Wtime() - runtime);
      printf("--------------------------------------------------!\n");
      printf("***************dt=%lf***************\n", dt);
      printf("***************flag=%d***************\n", flag1);
      printf("--------------------------------------------------!\n");
    }
    iter++;
    gyq++;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  walltime = MPI_Wtime() - walltime;
  MPI_Reduce(&walltime, &runtime, 1, MPI_DOUBLE, MPI_MAX, prank, MPI_COMM_WORLD);
  if (myrank == prank)
  {
    printf("time\t\t%lf\n", tt);
    printf("wall time\t%lf\n", runtime);
  }
  check_soln_new(tt);
}

void anisotropic_input(void)
{
  int i, j, n;
  double D[3][3];
  char filename[1024];
  FILE *fp;
  D[0][0] = kkx;
  D[0][1] = 0.0;
  D[0][2] = 0.0;
  D[1][0] = 0.0;
  D[1][1] = kky;
  D[1][2] = 0.0;
  D[2][0] = 0.0;
  D[2][1] = 0.0;
  D[2][2] = kkz;

  sprintf(filename, "anisotropic_input.txt");
  fp = fopen(filename, "r");
  for (n = 0; n < nac; n++)
  {
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        fscanf(fp, "%lf", &(ac[n].lambda[i][j]));
      }
    }
  }
  fclose(fp);

  for (n = 0; n < nac; n++)
  {
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        ac[n].lambda[i][j] = D[i][j] - ac[n].lambda[i][j];
      }
    }
  }
}

// size是指将数组压缩至char之后的大小
void float2bin(float *input, unsigned char *buff, int size)
{
  unsigned char temp_int = 1;
  unsigned char temp_c = 0;
  for (int temp_i = 0; temp_i < size; temp_i++)
  {
    temp_int = 1;
    temp_c = 0;
    for (int k = 0; k < 8; k++)
    {
      if (input[temp_i * 8 + k] > 0.1)
      {
        temp_c = temp_c | temp_int;
        temp_int = temp_int << 1;
        // temp_c = temp_c << 1;
        // printf("k:%d %d yes\n", k, temp_c);
      }
      else
      {
        temp_int = temp_int << 1;
        // printf("k:%d %d\n", k, temp_c);
      }
    }

    buff[temp_i] = temp_c;
  }
}

void compress_save_output(char *output_path, float *input, unsigned char *output, int size)
{
  float2bin(input, output, size);

  FILE *fp_out = fopen(output_path, "w"); //打开输出文件
  if (fp_out == NULL)
  {
    printf("Fail to save\n");
    exit(0);
  }

  for (int i = 0; i < size; i++)
  // for (int i = 0; i < 64*64*64*32; i++)
  {
    // printf("%f\n", (output_buffer + ((buffer_i + 1) % 2) * buffer_block_size)[i]);
    // fprintf(fp_out, "%f\n", (output_buffer + ((buffer_i + 1) % 2) * buffer_block_size)[i]);
    fprintf(fp_out, "%d\n", (output)[i]); // res_4_conv2
    // if ((output)[i] != 0)
    //   printf("%d\n", (output)[i]);
  }
  fclose(fp_out);
}
