#ifndef _SCLETD_H_
#define _SCLETD_H_

#include "mpi.h"

#define ScLETD_DEBUG
#define MODE 3
#define BLOCK 64

#ifdef ScLETD_DEBUG
#define ScLETDdebugTime MPI_Wtime()
#define ScLETDdebugPrintf printf
#else
#define ScLETDdebugTime 0
#define ScLETDdebugPrintf \
  while (0)               \
  printf
#endif

#define PI 3.141592653

int ioutput;
int ELASTIC;
int VariableMobility;
double ElasticScale;
int ANISOTROPIC;
int ETD2;
int *fieldgx, *fieldgy, *fieldgz;
int nac, nch;
int iter;
int left, right, top, bottom, front, back;
int nx, ny, nz;
int nprocs, procs[3], myrank, prank, cart_id[3];
int periodic;
//double PI;
double *MPX, *MPY, *MPZ, *MPXI, *MPYI, *MPZI, *DDX, *DDY, *DDZ;
double *MPX_b, *MPY_b, *MPZ_b, *MPXI_b, *MPYI_b, *MPZI_b, *DDX_b, *DDY_b, *DDZ_b;
double A1, A2, A3, A4, A5, C1;
double dt, t_total;
double epn2;
double xmin, ymin, zmin;
double xmax, ymax, zmax;
double hx, hy, hz;
double kkx, kky, kkz;
int restart;
char work_dir[1024];
char data_dir[1024];
int restart_iter;
char restart_dir[1024];
int ini_num;
int nghost;
int ix1, ix2, ix3, ix4, iy1, iy2, iy3, iy4, iz1, iz2, iz3, iz4;
int lnx, lny, lnz, gnx, gny, gnz;
double alpha, beta;
int stage;
MPI_Datatype left_right, top_bottom, front_back;
MPI_Status *status;
MPI_Comm YZ_COMM;

struct Allen_Cahn
{
  double u;
  double LE, KE;
  double lambda[3][3];
  double *f1, *f2, *f3;
  double *fieldE, *fieldE1, *fieldE2, *fieldEt, *fieldEp, *fieldE1p;
  double *fieldEs_left, *fieldEr_left, *fieldEs_right, *fieldEr_right;
  double *fieldEs_top, *fieldEr_top, *fieldEs_bottom, *fieldEr_bottom;
  double *fieldEr_front, *fieldEr_back;
  double *fieldEe_left, *fieldEe_right, *fieldEe_top, *fieldEe_bottom, *fieldEe_front, *fieldEe_back;
  double *fieldEu_left, *fieldEu_right, *fieldEu_top, *fieldEu_bottom, *fieldEu_front, *fieldEu_back;
  double *fieldEmu_left, *fieldEmu_right, *fieldEmu_top, *fieldEmu_bottom, *fieldEmu_front, *fieldEmu_back;
  double *phiE, *phiE2;
  double *felas;
  double *gradx, *grady, *gradz;
  MPI_Request *ireq_left_right_fieldE, *ireq_top_bottom_fieldE, *ireq_front_back_fieldE;
} * ac;

struct Cahn_Hilliard
{
  double c;
  double LCI, KCI, C;
  double *f0, *f1, *f2, *ft;
  double *M;
  double *fieldCI, *fieldCI1, *fieldCI2, *fieldCIt, *fieldCIp, *fieldCI1p;
  double *fieldCIs_left, *fieldCIr_left, *fieldCIs_right, *fieldCIr_right;
  double *fieldCIs_top, *fieldCIr_top, *fieldCIs_bottom, *fieldCIr_bottom;
  double *fieldCIr_front, *fieldCIr_back;
  double *fieldCIe_left, *fieldCIe_right, *fieldCIe_top, *fieldCIe_bottom, *fieldCIe_front, *fieldCIe_back;
  double *fieldCIu_left, *fieldCIu_right, *fieldCIu_top, *fieldCIu_bottom, *fieldCIu_front, *fieldCIu_back;
  double *fieldCImu_left, *fieldCImu_right, *fieldCImu_top, *fieldCImu_bottom, *fieldCImu_front, *fieldCImu_back;
  double *phiCI, *phiCI2;
  double *felas;
  double *felase_left, *felase_right, *felase_top, *felase_bottom, *felase_front, *felase_back;
  MPI_Request *ireq_left_right_fieldCI, *ireq_top_bottom_fieldCI, *ireq_front_back_fieldCI;
} * ch;

enum TYPE
{
  NEITHER_FUNCTION, AC_FUNCTION, CH_FUNCTION, BOTH_FUNCTION
};
 
void init_KL();
double EF();
double AC_f(int n);
double CH_f(int n);
// L(-F(C,E)+KC)
#define SWTICH_AC_FIELD1(f)                              \
  {                                                      \
    for (m = 0; m < nac; m++)                            \
    {                                                    \
      ac[m].u = ac[m].fieldE[k * nx * ny + j * nx + i];  \
    }                                                    \
    for (m = 0; m < nch; m++)                            \
    {                                                    \
      ch[m].c = ch[m].fieldCI[k * nx * ny + j * nx + i]; \
    }                                                    \
    f = ac[n].LE * ac[n].KE * ac[n].u;                   \
    f -= ac[n].LE * AC_f(n);                             \
  }

// L(F(C,E)-KC)
#define SWITCH_CH_FIELD1(f)                              \
  {                                                      \
    for (m = 0; m < nac; m++)                            \
    {                                                    \
      ac[m].u = ac[m].fieldE[k * nx * ny + j * nx + i];  \
    }                                                    \
    for (m = 0; m < nch; m++)                            \
    {                                                    \
      ch[m].c = ch[m].fieldCI[k * nx * ny + j * nx + i]; \
    }                                                    \
    f = -ch[n].LCI * ch[n].KCI * ch[n].c;                \
    f += ch[n].LCI * CH_f(n);                            \
  }

#define NAME(prefix, suffix) prefix##suffix
// ? - L(F(C,E))
#define SWITCH_CH_FIELDMU(suffix)                               \
  {                                                             \
    for (m = 0; m < nac; m++)                                   \
    {                                                           \
      ac[m].u = ac[m].NAME(fieldEe_, suffix)[l_e];              \
    }                                                           \
    for (m = 0; m < nch; m++)                                   \
    {                                                           \
      ch[m].c = ch[m].NAME(fieldCIe_, suffix)[l_e];             \
    }                                                           \
    NAME(fieldmu_, suffix)                                      \
    [l_mu] = ch[n].LCI * epn2 * NAME(fieldmu_, suffix)[l_mu];   \
    NAME(fieldmu_, suffix)                                      \
    [l_mu] -= ch[n].LCI * (CH_f(n));                            \
    if (ELASTIC == CH_FUNCTION)                                 \
    {                                                           \
      NAME(fieldmu_, suffix)                                    \
      [l_mu] -= ch[n].LCI * (ch[n].NAME(felase_, suffix)[l_e]); \
    }                                                           \
  }
#endif
