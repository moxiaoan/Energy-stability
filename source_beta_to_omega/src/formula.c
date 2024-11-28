#include <stdio.h>
#include "mpi.h"
#include "ScLETD.h"
#include <string.h>
#include <errno.h>
#include <stdlib.h>


void init_KL()
{
  ac[0].LE = 1.0;
  ac[1].LE = 1.0;
  ac[2].LE = 1.0;
  ac[3].LE = 1.0;
  ac[0].KE = 42.0;
  ac[1].KE = 42.0;
  ac[2].KE = 42.0;
  ac[3].KE = 42.0;
}

double
EF()
{
  return 0;
}

#define f_A (2376.0 / 2372.0) // 1
#define f_B ((7128.0 + 12.0 * 2372.0) / 2372.0) // 15
#define f_C ((4752.0 + 12.0 * 2372.0) / 2372.0) // 14
inline double AC_f(int n)
{
  double f;
  double sumh22 = ac[0].u * ac[0].u + ac[1].u * ac[1].u + ac[2].u * ac[2].u + ac[3].u * ac[3].u;
  f = f_A * ac[n].u - f_B * ac[n].u * ac[n].u + f_C * sumh22 * ac[n].u;
  //grain boundary
  //f += 66.0 * ac[n].u; 
  return f;
}

inline double CH_f(int n)
{
  double f;
  switch (n)
  {
  case 0:
  {
    break;
  }
  default:
  {
    printf("error:please choose correct CH_FIELD function.\n");
    exit(1);
  }
  }
  return f;
}
