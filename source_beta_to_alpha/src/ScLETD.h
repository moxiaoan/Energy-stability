#ifndef _SCLETD_H_
#define _SCLETD_H_

#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <mm_malloc.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include <complex.h>
#include "add.h"


#define LEN 256
#define MODE 3
#define BLOCK 64
#define THREADS_PER_BLOCK_X 4
#define THREADS_PER_BLOCK_Y 4
#define THREADS_PER_BLOCK_Z 4

#define PI 3.141592653589793
#define T_a 1093.0
//#define T 1173.0
#define R_a 8.31451
#define gnormal 50000
#define Bnormal 1.0e-18
#define ALPHA 0.02 // 1.07e2

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))

typedef double Dtype;
typedef unsigned short int Stype;
typedef unsigned char uint8_t;

#define ORI 7
#define DIM 3
// elastic
#define NELASTIC 6
int P, Q, count_gyq;
double Gbcc, Ghcp;
int ot2, ot1;
double *theta;
int flag, flag1;
int Approx;
int NX, NY, NZ;
double epsilon2d[NELASTIC][DIM * DIM];
double sigma2d[NELASTIC][DIM * DIM];
double C4d[DIM * DIM * DIM * DIM];
double *Bn[NELASTIC][NELASTIC], *Bn1[NELASTIC][NELASTIC];
double *tmpy_re, *tmpy_im, *tmpy_fftre[NELASTIC], *tmpy_fftim[NELASTIC], *Bn_fftre[NELASTIC], *Bn_fftim[NELASTIC], *fftre[NELASTIC], *fftim[NELASTIC];
double *s_left, *s_right, *s_top, *s_bottom, *s_front, *s_back;
double *r_left, *r_right, *r_top, *r_bottom, *r_front, *r_back;
double *e_left, *e_right, *e_top, *e_bottom, *e_front, *e_back;
MPI_Request *ireq_left_right, *ireq_top_bottom, *ireq_front_back;

double *Theta;
double *Epsilon2d[NELASTIC];
double *Sigma2d[NELASTIC];
double *omega_cpu;
double *C4D;
double *n11, *n22, *n33, *n0;
double *n123;
double *tmp0[NELASTIC][NELASTIC], *tmp1[NELASTIC][NELASTIC], *tmp2[NELASTIC][NELASTIC];
double *omega_inver;
double *E;
double *omega;
double *BN[NELASTIC][NELASTIC], *BN1[NELASTIC][NELASTIC];
double *tmpy_RE1, *tmpy_fftRE[NELASTIC], *tmpy_fftIM[NELASTIC], *BN_fftRE[NELASTIC], *BN_fftIM[NELASTIC], *fftRE[NELASTIC], *fftIM[NELASTIC];
double *Elas, *elas;
double *S_left, *S_right, *S_top, *S_bottom, *S_front, *S_back;
double *R_left, *R_right, *R_top, *R_bottom, *R_front, *R_back;

float *field2_all;
unsigned char *field2_cp;

double ElasticScale;
int ELASTIC;
int ANISOTROPIC;
int checkpoint_; // 不能叫checkpoint, 否则会出现重复定义(但是重复的地方不知道在哪也许是tf库里)

Dtype wmega, GHSERAL, GHSERTI, GHSERV, GBCCAL, GBCCTI, GHCPAL, GHCPV,
    BL13_0, BL13_1, BL13_2, BL132_0, BL12_0, BL12_1, BL32_0, BL32_1, BL32_2,
    HL13_0, HL13_1, HL13_2, HL132_0, HL132_1, HL132_2, HL12_0, HL12_1, HL32_0;
Dtype QALTI1, QALTI2, FALTI1, FALTI2, QALAL1, FALAL1, QALV1, FALV1, QTITI1, QTITI2,
    FTITI1, FTITI2, QVTI1, QVTI2, FVTI1, FVTI2, QTIV1, QTIV2, FTIV1, FTIV2,
    QVV1, QVV2, FVV1, FVV2, ALTITIV0, ALTITIV1, ALVTIV0, ALVTIV1, ALTIALTI0,
    ALTIALTI1, ALALALTI0, ALALALTI1, ALALALV0, ALVALV0, ALALTIV0, VVAL12, VVAL14,
    DTALTI, DTALAL, ATALV, VV12, VV14, DTTITI, VV22, VV24, DTTIV, VV32, VV34, DTVTI, VV42, VV44, DTVV,
    AMTITIHCP, AMALTIHCP, AMVTIHCP, DTALV;
Dtype G1_1, G1_2, G1_3, G1_13_0, G1_13_1, G1_32_0, G1_12_0,
    G2_1, G2_2, G2_3, G2_12_0, G2_32_0, G2_32_1,
    G3_1, G3_2, G3_3, G3_13_0, G3_13_1, G3_32_0, G3_32_1,
    HG1_1, HG1_2, HG1_3, HG2_1, HG2_2, HG2_3, HG3_1, HG3_2, HG3_3;

int ioutput;
Dtype chcp2, cbcc2;
int restart, restart_iter;
int nout;
int nac, nch;
int iter;
int tt, out_tt, out_iter;
int left, right, top, bottom, front, back;
int nx, ny, nz;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int namelen;
int deviceId;

int nprocs, procs[3], myrank, prank, cart_id[3];
int periodic;
Dtype *MPX, *MPY, *MPZ, *MPXI, *MPYI, *MPZI, *DDX, *DDY, *DDZ;
Dtype *MPX_b, *MPY_b, *MPZ_b, *MPXI_b, *MPYI_b, *MPZI_b, *DDX_b, *DDY_b, *DDZ_b;
Dtype dt, t_total;
Dtype xmin, ymin, zmin;
Dtype xmax, ymax, zmax;
// new
int lr_size, tb_size, fb_size;

Dtype hx, hy, hz;
Dtype kkx, kky, kkz;
char work_dir[1024];
char data_dir[1024];
int ini_num;
int nghost;
int ix1, ix2, ix3, ix4, iy1, iy2, iy3, iy4, iz1, iz2, iz3, iz4;
int lnx, lny, lnz, gnx, gny, gnz;
Dtype alpha, beta, alpha1, beta1, beta2;
MPI_Datatype left_right, top_bottom, front_back;
MPI_Status *status;
MPI_Comm XYZ_COMM;
MPI_Comm YZ_COMM;
int irun, nchk, chk;
// offset
size_t offset, offset_Er;
size_t u_fb, u_tb, u_lr;
size_t e_fb, e_tb, e_lr;
// elastic offset
size_t ela_size;
// time

float dgemm_time, F1_time, F2_time, updateU_new_time, calc_time;
float trans_unpack_time;
float trans_mu_time;
float trans_enlarge_time;
float trans_pack_time;
float trans_Memcpy_time;
float trans_MPI_time;
float trans_time;
float timer;
float timer2;

// cpu variables
struct Allen_Cahn
{
  Dtype LE, KE;
  Dtype u;
  Dtype lambda[3][3];
  Dtype epn2;
  Stype *field2b;
  Dtype *fieldE, *fieldE_old;
  Dtype *fieldEs_left, *fieldEr_left, *fieldEs_right, *fieldEr_right;
  Dtype *fieldEs_top, *fieldEr_top, *fieldEs_bottom, *fieldEr_bottom;
  Dtype *fieldEr_front, *fieldEr_back;
  Dtype *fieldEs_front, *fieldEs_back;
  Dtype *fieldEe_left, *fieldEe_right, *fieldEe_top, *fieldEe_bottom, *fieldEe_front, *fieldEe_back;
  MPI_Request *ireq_left_right_fieldE, *ireq_top_bottom_fieldE, *ireq_front_back_fieldE;
} * ac;
double *lambdar, *lambdar1;

struct Cahn_Hilliard
{
  Dtype LE, KE;
  Dtype u;
  Dtype epn2;
  Stype *field2b;
  Dtype *fieldCI, *fieldCI_old;
  Dtype *fieldCIs_left, *fieldCIr_left, *fieldCIs_right, *fieldCIr_right;
  Dtype *fieldCIs_top, *fieldCIr_top, *fieldCIs_bottom, *fieldCIr_bottom;
  Dtype *fieldCIr_front, *fieldCIr_back;
  Dtype *fieldCIs_front, *fieldCIs_back;
  Dtype *fieldCIe_left, *fieldCIe_right, *fieldCIe_top, *fieldCIe_bottom, *fieldCIe_front, *fieldCIe_back;
  Dtype *fieldCImu_left, *fieldCImu_right, *fieldCImu_top, *fieldCImu_bottom, *fieldCImu_front, *fieldCImu_back;
  MPI_Request *ireq_left_right_fieldCI, *ireq_top_bottom_fieldCI, *ireq_front_back_fieldCI;
} * ch;

Stype *field2BE, *field2BCI;
Dtype *lambda, *f_aniso;
Dtype *fieldE, *fieldE_old, *fieldCI, *fieldCI_old, *fieldEr;
Dtype *fieldEu_top, *fieldEu_bottom, *fieldEu_front, *fieldEu_back, *fieldEu_left, *fieldEu_right;
Dtype *fieldCImu_top, *fieldCImu_bottom, *fieldCImu_front, *fieldCImu_back, *fieldCImu_left, *fieldCImu_right;
Dtype *fieldEr_front, *fieldCIr_front;
Dtype *fieldEr_back, *fieldCIr_back;
Dtype *fieldEr_top, *fieldEs_top, *fieldCIr_top, *fieldCIs_top;
Dtype *fieldEr_bottom, *fieldEs_bottom, *fieldCIr_bottom, *fieldCIs_bottom;
Dtype *fieldEr_left, *fieldEs_left, *fieldCIr_left, *fieldCIs_left;
Dtype *fieldEr_right, *fieldEs_right, *fieldCIr_right, *fieldCIs_right;
Dtype *fieldEe_front, *fieldEe_back, *fieldCIe_front, *fieldCIe_back;
Dtype *fieldEe_top, *fieldEe_bottom, *fieldCIe_top, *fieldCIe_bottom;
Dtype *fieldEe_left, *fieldEe_right, *fieldCIe_left, *fieldCIe_right;

Dtype *ddx, *ddy, *ddz;
Dtype *fieldE1, *fieldCI1;
Dtype *fieldEt, *fieldEp, *fieldCIt, *fieldCIp;
Dtype *dfc1, *dfc2, *ft, *f1, *f2, *f3;
Dtype *phiE, *phiCI;
Dtype *M, *C;
Dtype *mpxi, *mpyi, *mpzi, *mpx, *mpy, *mpz;

Dtype *Gx_re, *Gx_im, *Gy_re, *Gy_im, *Gz_re, *Gz_im, *temp1, *temp2, *temp3, *temp4, *temp5, *temp6, *temp7, *temp8, *Itemp1, *Itemp2, *Itemp3, *Itemp4, *Itemp5;

int *ori, *Ori, *f, *F;
double *field1, *field2;

double CH_f(int n);
void write_chk(void);
void write_field2B(int irun, int mode);
void read_chk(void);

#define NAME(prefix, suffix) prefix##suffix
// epn2 * laplace(ch.u) - f(ac.u, ch.u)
#define SWITCH_CH_UMU(suffix)                                  \
  {                                                            \
    int m;                                                     \
    for (m = 0; m < nac; m++)                                  \
    {                                                          \
      ac[m].u = ac[m].NAME(fieldEe_, suffix)[l_e];             \
    }                                                          \
    for (m = 0; m < nch; m++)                                  \
    {                                                          \
      ch[m].u = ch[m].NAME(fieldCIe_, suffix)[l_e];            \
    }                                                          \
    NAME(Umu_, suffix)                                         \
    [l_mu] = ch[n].LE * ch[n].epn2 * NAME(Umu_, suffix)[l_mu]; \
    NAME(Umu_, suffix)                                         \
    [l_mu] -= ch[n].LE * CH_f(n);                              \
  }

#endif
