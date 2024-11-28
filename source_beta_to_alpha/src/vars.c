#include "ScLETD.h"
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include "anisotropic_hip.h"

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
    gny = lny * procs[0] - 2 * nghost * procs[0];
  }
  else
  {
    gny = lny * procs[0] - 2 * nghost * (procs[0] - 1);
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
    gnz = lnz * procs[0] - 2 * nghost * procs[0];
  }
  else
  {
    gnz = lnz * procs[0] - 2 * nghost * (procs[0] - 1);
  }

  hx = (xmax - xmin) / gnx;
  hy = (ymax - ymin) / gny;
  hz = (zmax - zmin) / gnz;
}

void init_vars(void)
{
  int i, j, k, l;

  alpha = 1.0;
  beta = 0.0;
  alpha1 = -1.0;
  beta1 = -1.0;
  beta2 = 1.0;

  ac[0].LE = 2.0;
  ac[1].LE = 2.0;
  ac[2].LE = 2.0;
  ac[3].LE = 2.0;
  ac[4].LE = 2.0;
  ac[5].LE = 2.0;
  ac[0].KE = 3.0;//30.0;//80.0;//175.0; // 175.0;//616.452;//616425.0;//175.0;//45.0;//280;//200;//175.0;//4.0;
  ac[1].KE = 3.0;//30.0;//80.0;//175.0; // 175.0;//616.452;//616425.0;//175.0;//45.0;//280;//200;//175.0;//4.0;
  ac[2].KE = 3.0;//30.0;//80.0;//175.0;
  ac[3].KE = 3.0;//30.0;//80.0;//175.0;
  ac[4].KE = 3.0;//30.0;//80.0;//175.0;
  ac[5].KE = 3.0;//30.0;//80.0;//175.0;
// ch[0].LE = 1.0;
// ch[1].LE = 1.0;
// ch[0].KE = 8.0; // 8.0;
// ch[1].KE = 8.0; // 8.0;

  ac[0].epn2 = 0.048;
  ac[1].epn2 = 0.048;
  ac[2].epn2 = 0.048;
  ac[3].epn2 = 0.048;
  ac[4].epn2 = 0.048;
  ac[5].epn2 = 0.048;
// ch[0].epn2 = 0.0;
// ch[1].epn2 = 0.0;

  wmega = 100.0;

  GHSERAL = -11277.683 + 188.661987 * T_a - 31.748192 * T_a * log(T_a) - 1.234E+28 * pow(T_a, (-9)) + 137529.83; // 8.9664e+04
  GHSERTI = +908.837 + 67.048538 * T_a - 14.9466 * T_a * log(T_a) - 0.0081465 * pow(T_a, 2) + 2.02715E-07 * pow(T_a, 3) - 1477660 * pow(T_a, (-1)) + 58166.677;
  // 7.2379e+03
  GHSERV = -7967.84 + 143.291 * T_a - 25.9 * T_a * log(T_a) + 6.25E-05 * pow(T_a, 2) - 6.8E-07 * pow(T_a, 3) + 77277.217; // 2.7046e+04

  GBCCAL = 10083 - 4.812 * T_a + GHSERAL; // 9.4487e+04
  GBCCTI = +5758.6 + 38.38987 * T_a - 7.4305 * T_a * log(T_a) + 0.0093636 * pow(T_a, 2) - 1.04805E-06 * pow(T_a, 3) - 525093 * pow(T_a, (-1)) + GHSERTI;
  // 7.4701e+03

  GHCPAL = +5481 - 1.799 * T_a + GHSERAL; // 9.3179e+04
  GHCPV = +4000 + 2.4 * T_a + GHSERV;     // 3.3669e+04

  // 3=TI 1=AL 2=V

  BL13_0 = -128500 + 39 * T_a; //-85873
  BL13_1 = 6000;
  BL13_2 = 21200;

  BL132_0 = -8000;

  BL12_0 = -98000 + 32 * T_a; //-63024
  BL12_1 = +3500 - 5 * T_a;   //-1965

  BL32_0 = +10500.0 - 1.5 * T_a; // 8860.5
  BL32_1 = 2000.0;
  BL32_2 = 1000.0;

  HL13_0 = -133750 + 39 * T_a; //-91123
  HL13_1 = 250;
  HL13_2 = 17250;

  HL132_0 = -60000;
  HL132_1 = -30000;
  HL132_2 = 60000;

  HL12_0 = -98000 + 32 * T_a; //-63024
  HL12_1 = +3500 - 5 * T_a;   //-1965

  HL32_0 = +30250 - 10.0 * T_a; // 19320

  QALTI1 = 204000;
  QALTI2 = 96000;
  FALTI1 = 5.51E-6;
  FALTI2 = 5.19E-10;

  QALAL1 = 19000;
  FALAL1 = 5E-4;

  QALV1 = 290000;
  FALV1 = 1.85E-5;

  QTITI1 = 121000;
  QTITI2 = 237000;
  FTITI1 = 1.47E-8;
  FTITI2 = 5.91E-5;

  QVTI1 = 107000;
  QVTI2 = 212000;
  FVTI1 = 1.55E-9;
  FVTI2 = 2.54E-5;

  QTIV1 = 320000;
  QTIV2 = 544000;
  FTIV1 = 1.62E-4;
  FTIV2 = 51;

  QVV1 = 306000;
  QVV2 = 493000;
  FVV1 = 2.73E-5;
  FVV2 = 1.58;

  ALTITIV0 = 111000 - 20 * T_a; // 89140
  ALTITIV1 = -43000 + 10 * T_a;
  ALVTIV0 = 78000 + 11 * T_a;
  ALVTIV1 = 47000 - 22 * T_a;

  ALTIALTI0 = -247000 - 83 * T_a;
  ALTIALTI1 = -146000;
  ALALALTI0 = -1059000 + 373 * T_a;
  ALALALTI1 = -921000 + 342 * T_a;

  ALALALV0 = -4143000;
  ALVALV0 = -5752000;

  ALALTIV0 = +180000 - 277 * T_a;

  VVAL12 = -0.12027 * QALTI1 * pow(T_a, (-1));
  VVAL14 = -0.12027 * QALTI2 * pow(T_a, (-1));
  DTALTI = FALTI1 * exp(VVAL12) + FALTI2 * exp(VVAL14);

  DTALAL = -QALAL1 + R_a * T_a * log(FALAL1);

  DTALV = -QALV1 + R_a * T_a * log(FALV1);

  VV12 = -0.12027 * QTITI1 * pow(T_a, -1);
  VV14 = -0.12027 * QTITI2 * pow(T_a, -1);
  DTTITI = FTITI1 * exp(VV12) + FTITI2 * exp(VV14);
  VV22 = -0.12027 * QTIV1 * pow(T_a, -1);
  VV24 = -0.12027 * QTIV2 * pow(T_a, -1);
  DTTIV = FTIV1 * exp(VV22) + FTIV2 * exp(VV24);
  VV32 = -0.12027 * QVTI1 * pow(T_a, -1);
  VV34 = -0.12027 * QVTI2 * pow(T_a, -1);
  DTVTI = FVTI1 * exp(VV32) + FVTI2 * exp(VV34);
  VV42 = -0.12027 * QVV1 * pow(T_a, -1);
  VV44 = -0.12027 * QVV2 * pow(T_a, -1);
  DTVV = FVV1 * exp(VV42) + FVV2 * exp(VV44);

  AMTITIHCP = -3.03E5 + R_a * T_a * log(1.35E-3);
  AMALTIHCP = -3.29E5 + R_a * T_a * log(6.6E-3);
  AMVTIHCP = -2.50E5 + R_a * T_a * log(1E-3);

  // mobility of Al (B2)
  G1_1 = DTALAL;
  G1_2 = DTALV;
  G1_3 = R_a * T_a * log(DTALTI);

  G1_13_0 = ALALALTI0;
  G1_13_1 = ALALALTI1;
  G1_32_0 = ALALTIV0;
  G1_12_0 = ALALALV0;

  // mobility of V (B2)
  G2_1 = DTALAL;
  G2_2 = R_a * T_a * log(DTVV);
  G2_3 = R_a * T_a * log(DTVTI);

  G2_12_0 = ALVALV0;
  G2_32_0 = ALVTIV0;
  G2_32_1 = ALVTIV1;

  // mobility of TI (B3)
  G3_1 = DTALAL;
  G3_2 = R_a * T_a * log(DTTIV);
  G3_3 = R_a * T_a * log(DTTITI);

  G3_13_0 = ALTIALTI0;
  G3_13_1 = ALTIALTI1;
  G3_32_0 = ALTITIV0;
  G3_32_1 = ALTITIV1;

  // mobility of Al (B2)
  HG1_1 = AMALTIHCP;
  HG1_2 = AMALTIHCP;
  HG1_3 = AMALTIHCP;

  // mobility of V (B2)
  HG2_1 = AMVTIHCP;
  HG2_2 = AMVTIHCP;
  HG2_3 = AMVTIHCP;

  // mobility of TI (B3)
  HG3_1 = AMTITIHCP;
  HG3_2 = AMTITIHCP;
  HG3_3 = AMTITIHCP;

  offset = nx * ny * nz;
  offset_Er = (nx + 2 * 2) * (ny + 2 * 2) * (nz + 2 * 2);
  lr_size = (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2);
  tb_size = nx * (nghost + 2) * (nz + (nghost + 2) * 2);
  fb_size = nx * ny * (nghost + 2);
  u_fb = nx * ny;
  u_tb = nx * nz;
  u_lr = ny * nz;
  e_fb = (nx + 4) * (ny + 4) * (nghost + 2);
  e_tb = (nx + 4) * (nghost + 2) * (nz + 4);
  e_lr = (nghost + 2) * (ny + 4) * (nz + 4);
  ela_size = Approx * Approx * Approx;
  x_m = nx;
  x_n = ny * nz;
  x_k = nx;
  y_m = ny;
  y_n = nz * nx;
  y_k = ny;
  z_m = nz;
  z_n = nx * ny;
  z_k = nz;

  // Gx_m = 16;
  Gx_m = nx - 2 * nghost;
  Gx_n = (ny - 2 * nghost) * (nz - 2 * nghost);
  Gx_k = Approx;
  Gy_m = ny - 2 * nghost;
  Gy_n = (nz - 2 * nghost) * Approx;
  Gy_k = Approx;
  Gz_m = nz - 2 * nghost;
  Gz_n = Approx * Approx;
  Gz_k = Approx;
}

void alloc_vars(void)
{
  int n;
  ac = (struct Allen_Cahn *)_mm_malloc(nac * sizeof(struct Allen_Cahn), 256);
  for (n = 0; n < nac; n++)
  {
    ac[n].field2b = (Stype *)_mm_malloc(sizeof(Stype) * nx * ny * nz, 256);
    ac[n].fieldE = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ac[n].fieldE_old = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
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
    ac[n].fieldEs_front = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ac[n].fieldEs_back = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ac[n].fieldEe_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ac[n].fieldEe_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ac[n].fieldEe_top = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ac[n].fieldEe_bottom = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ac[n].fieldEe_front = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ac[n].fieldEe_back = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ac[n].ireq_left_right_fieldE = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ac[n].ireq_top_bottom_fieldE = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ac[n].ireq_front_back_fieldE = (MPI_Request *)calloc(4, sizeof(MPI_Request));
  }
  ori = (int *)_mm_malloc(sizeof(int) * nx * ny * nz, 256);
  f = (int *)_mm_malloc(sizeof(int) * nx * ny * nz, 256);
  field2 = (double *)_mm_malloc(sizeof(double) * (nx - 2 * nghost) * (ny - 2 * nghost) * (nz - 2 * nghost) / 8, 256);

  field2_all = (float *)_mm_malloc(sizeof(float) * 128 * 128 * 128 * nac, 256);
  field2_cp = (unsigned char *)_mm_malloc(sizeof(unsigned char) * 128 * 128 * 16, 256);
  ch = (struct Cahn_Hilliard *)_mm_malloc(nch * sizeof(struct Cahn_Hilliard), 256);
  for (n = 0; n < nch; n++)
  {
    ch[n].field2b = (Stype *)_mm_malloc(sizeof(Stype) * nx * ny * nz, 256);
    ch[n].fieldCI = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
    ch[n].fieldCI_old = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
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
    ch[n].fieldCIs_front = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ch[n].fieldCIs_back = (double *)_mm_malloc(sizeof(double) * nx * ny * (nghost + 2), 256);
    ch[n].fieldCIe_left = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ch[n].fieldCIe_right = (double *)_mm_malloc(sizeof(double) * (nghost + 2) * (ny + 4) * (nz + 4), 256);
    ch[n].fieldCIe_top = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ch[n].fieldCIe_bottom = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (nghost + 2) * (nz + 4), 256);
    ch[n].fieldCIe_front = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ch[n].fieldCIe_back = (double *)_mm_malloc(sizeof(double) * (nx + 4) * (ny + 4) * (nghost + 2), 256);
    ch[n].fieldCImu_left = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ch[n].fieldCImu_right = (double *)_mm_malloc(sizeof(double) * ny * nz, 256);
    ch[n].fieldCImu_top = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ch[n].fieldCImu_bottom = (double *)_mm_malloc(sizeof(double) * nx * nz, 256);
    ch[n].fieldCImu_front = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);
    ch[n].fieldCImu_back = (double *)_mm_malloc(sizeof(double) * nx * ny, 256);

    ch[n].ireq_left_right_fieldCI = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ch[n].ireq_top_bottom_fieldCI = (MPI_Request *)calloc(4, sizeof(MPI_Request));
    ch[n].ireq_front_back_fieldCI = (MPI_Request *)calloc(4, sizeof(MPI_Request));
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
  rocblas_create_handle(&handle);

  hipMalloc(&field1, ((nx - 2 * nghost) * (ny - 2 * nghost) * (nz - 2 * nghost) / 8) * sizeof(Dtype));
  hipMalloc(&f1, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&f2, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&f3, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&ft, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&dfc1, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&dfc2, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldE, nac * nx * ny * nz * sizeof(Dtype));
  // hipMalloc (&fieldE_old, nac * nx * ny * nz * sizeof (Dtype));
  hipMalloc(&fieldE1, nac * nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldEt, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldEp, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&phiE, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldCI, nch * nx * ny * nz * sizeof(Dtype));
  // hipMalloc (&fieldCI_old, nch * nx * ny * nz * sizeof (Dtype));
  hipMalloc(&fieldCI1, nch * nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldCIt, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldCIp, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&phiCI, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&M, 2 * nch * nx * ny * nz * sizeof(Dtype));
  hipMalloc(&C, nch * sizeof(Dtype));
  hipMalloc (&field2BE, nac * nx * ny * nz * sizeof (Stype));
  hipMalloc (&field2BCI, nch * nx * ny * nz * sizeof (Stype));

  hipMalloc(&fieldEu_left, nac * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldEu_right, nac * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldEu_front, nac * nx * ny * sizeof(Dtype));
  hipMalloc(&fieldEu_back, nac * nx * ny * sizeof(Dtype));
  hipMalloc(&fieldEu_top, nac * nx * nz * sizeof(Dtype));
  hipMalloc(&fieldEu_bottom, nac * nx * nz * sizeof(Dtype));

  hipMalloc(&fieldCImu_left, nch * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldCImu_right, nch * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldCImu_front, nch * nx * ny * sizeof(Dtype));
  hipMalloc(&fieldCImu_back, nch * nx * ny * sizeof(Dtype));
  hipMalloc(&fieldCImu_top, nch * nx * nz * sizeof(Dtype));
  hipMalloc(&fieldCImu_bottom, nch * nx * nz * sizeof(Dtype));
  hipMalloc(&lambda, nac * 3 * 3 * sizeof(Dtype));
  hipMalloc(&f_aniso, nx * ny * nz * sizeof(Dtype));
  hipMalloc(&fieldEr, nac * (nx + 2 * 2) * (ny + 2 * 2) * (nz + 2 * 2) * sizeof(Dtype));
  hipMalloc(&Ori, nx * ny * nz * sizeof(int));
  hipMalloc(&F, nx * ny * nz * sizeof(int));

  hipMalloc(&mpxi, nx * nx * sizeof(Dtype));
  hipMalloc(&mpyi, ny * ny * sizeof(Dtype));
  hipMalloc(&mpzi, nz * nz * sizeof(Dtype));
  hipMalloc(&mpx, nx * nx * sizeof(Dtype));
  hipMalloc(&mpy, ny * ny * sizeof(Dtype));
  hipMalloc(&mpz, nz * nz * sizeof(Dtype));
  hipMalloc(&ddx, nx * sizeof(Dtype));
  hipMalloc(&ddy, ny * sizeof(Dtype));
  hipMalloc(&ddz, nz * sizeof(Dtype));

  hipMalloc(&fieldEr_left, nac * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEs_left, nac * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEr_right, nac * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEs_right, nac * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEr_top, nac * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEs_top, nac * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEr_bottom, nac * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEs_bottom, nac * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldEr_front, nac * nx * ny * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldEr_back, nac * nx * ny * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldEe_front, nac * (nx + 4) * (ny + 4) * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldEe_back, nac * (nx + 4) * (ny + 4) * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldEe_top, nac * (nx + 4) * (nghost + 2) * (nz + 4) * sizeof(Dtype));
  hipMalloc(&fieldEe_bottom, nac * (nx + 4) * (nghost + 2) * (nz + 4) * sizeof(Dtype));
  hipMalloc(&fieldEe_left, nac * (nghost + 2) * (ny + 4) * (nz + 4) * sizeof(Dtype));
  hipMalloc(&fieldEe_right, nac * (nghost + 2) * (ny + 4) * (nz + 4) * sizeof(Dtype));

  hipMalloc(&fieldCIr_left, nch * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIs_left, nch * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIr_right, nch * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIs_right, nch * (nghost + 2) * (ny + (nghost + 2) * 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIr_top, nch * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIs_top, nch * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIr_bottom, nch * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIs_bottom, nch * nx * (nghost + 2) * (nz + (nghost + 2) * 2) * sizeof(Dtype));
  hipMalloc(&fieldCIr_front, nch * nx * ny * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldCIr_back, nch * nx * ny * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldCIe_front, nch * (nx + 4) * (ny + 4) * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldCIe_back, nch * (nx + 4) * (ny + 4) * (nghost + 2) * sizeof(Dtype));
  hipMalloc(&fieldCIe_top, nch * (nx + 4) * (nghost + 2) * (nz + 4) * sizeof(Dtype));
  hipMalloc(&fieldCIe_bottom, nch * (nx + 4) * (nghost + 2) * (nz + 4) * sizeof(Dtype));
  hipMalloc(&fieldCIe_left, nch * (nghost + 2) * (ny + 4) * (nz + 4) * sizeof(Dtype));
  hipMalloc(&fieldCIe_right, nch * (nghost + 2) * (ny + 4) * (nz + 4) * sizeof(Dtype));

  hipEventCreate(&st);
  hipEventCreate(&ed);
  hipEventCreate(&st2);
  hipEventCreate(&ed2);
}
void read_matrices(void)
{
  int n;
  char filename[1024], dirname[1024], wirname[1024];
  int i, j, k, l, id, mp_ofst, d_ofst;
  FILE *file;

  if (cart_id[0] == 0 and cart_id[1] == 0)
  {
    sprintf(dirname, "%s%02d", data_dir, cart_id[2]);
    if (mkdir(dirname, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH) < 0)
    {
      // printf("mkdir failed\n");
      // return 2;
    }
  }
  if (cart_id[0] == 0 and cart_id[1] == 0)
  {
    sprintf(wirname, "%s%02d", work_dir, cart_id[2]);
    if (mkdir(wirname, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH) < 0)
    {
      // printf("mkdir failed\n");
      // return 2;
    }
  }

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
}

void dealloc_vars(void)
{
  int l;
  int n;
  for (n = 0; n < nac; n++)
  {
    _mm_free(ac[n].field2b);
    _mm_free(ac[n].fieldE);
    _mm_free(ac[n].fieldE_old);
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
    _mm_free(ac[n].fieldEs_front);
    _mm_free(ac[n].fieldEs_back);
    _mm_free(ac[n].fieldEe_left);
    _mm_free(ac[n].fieldEe_right);
    _mm_free(ac[n].fieldEe_top);
    _mm_free(ac[n].fieldEe_bottom);
    _mm_free(ac[n].fieldEe_front);
    _mm_free(ac[n].fieldEe_back);

    for (l = 0; l < 4; l++)
    {
      MPI_Request_free(&ac[n].ireq_left_right_fieldE[l]);
      MPI_Request_free(&ac[n].ireq_top_bottom_fieldE[l]);
      MPI_Request_free(&ac[n].ireq_front_back_fieldE[l]);
    }
    free(ac[n].ireq_left_right_fieldE);
    free(ac[n].ireq_top_bottom_fieldE);
    free(ac[n].ireq_front_back_fieldE);
  }
  _mm_free(ac);
  _mm_free(ori);
  _mm_free(f);
  _mm_free(field2);

  for (n = 0; n < nch; n++)
  {
    _mm_free(ch[n].field2b);
    _mm_free(ch[n].fieldCI);
    _mm_free(ch[n].fieldCI_old);
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
    _mm_free(ch[n].fieldCIs_front);
    _mm_free(ch[n].fieldCIs_back);
    _mm_free(ch[n].fieldCIe_left);
    _mm_free(ch[n].fieldCIe_right);
    _mm_free(ch[n].fieldCIe_top);
    _mm_free(ch[n].fieldCIe_bottom);
    _mm_free(ch[n].fieldCIe_front);
    _mm_free(ch[n].fieldCIe_back);
    _mm_free(ch[n].fieldCImu_left);
    _mm_free(ch[n].fieldCImu_right);
    _mm_free(ch[n].fieldCImu_top);
    _mm_free(ch[n].fieldCImu_bottom);
    _mm_free(ch[n].fieldCImu_front);
    _mm_free(ch[n].fieldCImu_back);

    for (l = 0; l < 4; l++)
    {
      MPI_Request_free(&ch[n].ireq_left_right_fieldCI[l]);
      MPI_Request_free(&ch[n].ireq_top_bottom_fieldCI[l]);
      MPI_Request_free(&ch[n].ireq_front_back_fieldCI[l]);
    }
    free(ch[n].ireq_left_right_fieldCI);
    free(ch[n].ireq_top_bottom_fieldCI);
    free(ch[n].ireq_front_back_fieldCI);
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
  //_mm_free(c);
  free(status);

  rocblas_destroy_handle(handle);
  hipFree(fieldEr_front);
  hipFree(fieldEr_back);
  hipFree(fieldEr_top);
  hipFree(fieldEs_top);
  hipFree(fieldEr_bottom);
  hipFree(fieldEs_bottom);
  hipFree(fieldEr_left);
  hipFree(fieldEs_left);
  hipFree(fieldEr_right);
  hipFree(fieldEs_right);

  hipFree(fieldCIr_front);
  hipFree(fieldCIr_back);
  hipFree(fieldCIr_top);
  hipFree(fieldCIs_top);
  hipFree(fieldCIr_bottom);
  hipFree(fieldCIs_bottom);
  hipFree(fieldCIr_left);
  hipFree(fieldCIs_left);
  hipFree(fieldCIr_right);
  hipFree(fieldCIs_right);

  hipFree(fieldEe_front);
  hipFree(fieldEe_back);
  hipFree(fieldEe_top);
  hipFree(fieldEe_bottom);
  hipFree(fieldEe_left);
  hipFree(fieldEe_right);

  hipFree(fieldEu_left);
  hipFree(fieldEu_right);
  hipFree(fieldEu_top);
  hipFree(fieldEu_bottom);
  hipFree(fieldEu_front);
  hipFree(fieldEu_back);

  hipFree(fieldCImu_left);
  hipFree(fieldCImu_right);
  hipFree(fieldCImu_top);
  hipFree(fieldCImu_bottom);
  hipFree(fieldCImu_front);
  hipFree(fieldCImu_back);

  hipFree(fieldCIe_front);
  hipFree(fieldCIe_back);
  hipFree(fieldCIe_top);
  hipFree(fieldCIe_bottom);
  hipFree(fieldCIe_left);
  hipFree(fieldCIe_right);

  hipFree(field1);
  hipFree(f1);
  hipFree(f2);
  hipFree(f3);
  hipFree(ft);
  hipFree(dfc1);
  hipFree(dfc2);
  hipFree(fieldE);
  // hipFree (fieldE_old);
  hipFree(fieldE1);
  hipFree(fieldEt);
  hipFree(fieldEp);
  hipFree(fieldCI);
  // hipFree (fieldCI_old);
  hipFree(fieldCI1);
  hipFree(fieldCIt);
  hipFree(fieldCIp);
  hipFree (field2BE);
  hipFree (field2BCI);
  hipFree(phiE);
  hipFree(M);
  hipFree(C);
  hipFree(phiCI);
  hipFree(f_aniso);
  hipFree(lambda);
  hipFree(fieldEr);
  hipFree(Ori);
  hipFree(F);

  hipFree(mpxi);
  hipFree(mpyi);
  hipFree(mpzi);
  hipFree(ddx);
  hipFree(ddy);
  hipFree(ddz);
  hipFree(mpx);
  hipFree(mpy);
  hipFree(mpz);

  hipEventDestroy(st);
  hipEventDestroy(ed);
  hipEventDestroy(st2);
  hipEventDestroy(ed2);
}

void write_chk(void)
{
  int i, j, k, n;
  char filename[1024];
  FILE *file;
  for (n = 0; n < nac; n++)
  {
    hipMemcpy(ac[n].fieldE, fieldE + n * nx * ny * nz, sizeof(Dtype) * nx * ny * nz, hipMemcpyDeviceToHost);

    sprintf(filename, "%s%02d/eta%d_chk%d_%02d%02d%02d.dat", work_dir, cart_id[2], n, chk, cart_id[0], cart_id[1], cart_id[2]);
    file = fopen(filename, "wb");

    fwrite(&irun, sizeof(int), 1, file);
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          fwrite(&ac[n].fieldE[k * nx * ny + j * nx + i], sizeof(Dtype), 1, file);
        }
      }
    }
    fclose(file);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    hipMemcpy(ch[n].fieldCI, fieldCI + n * nx * ny * nz, sizeof(Dtype) * nx * ny * nz, hipMemcpyDeviceToHost);

    sprintf(filename, "%s%02d/c%d_chk%d_%02d%02d%02d.dat", work_dir, cart_id[2], n, chk, cart_id[0], cart_id[1], cart_id[2]);
    file = fopen(filename, "wb");

    fwrite(&irun, sizeof(int), 1, file);
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          fwrite(&ch[n].fieldCI[k * nx * ny + j * nx + i], sizeof(Dtype), 1, file);
        }
      }
    }
    fclose(file);
  }
#endif
  // MPI_Barrier (MPI_COMM_WORLD);
}

void read_chk(void)
{
  int i, j, k, n;
  char filename[1024];
  FILE *file;

  for (n = 0; n < nac; n++)
  {
    sprintf(filename, "%s%02d/eta%d_chk%d_%02d%02d%02d.dat", work_dir, cart_id[2], n, chk, cart_id[0], cart_id[1], cart_id[2]);
    file = fopen(filename, "rb");

    fread(&irun, sizeof(int), 1, file);
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          fread(&ac[n].fieldE[k * nx * ny + j * nx + i], sizeof(Dtype), 1, file);
        }
      }
    }
    fclose(file);

    hipMemcpy(fieldE + n * nx * ny * nz, ac[n].fieldE, sizeof(Dtype) * nx * ny * nz, hipMemcpyHostToDevice);
  }
#if 0
  for (n = 0; n < nch; n++)
  {
    sprintf(filename, "%s%02d/c%d_chk%d_%02d%02d%02d.dat", work_dir, cart_id[2], n, chk, cart_id[0], cart_id[1], cart_id[2]);
    file = fopen(filename, "rb");

    fread(&irun, sizeof(int), 1, file);
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          fread(&ch[n].fieldCI[k * nx * ny + j * nx + i], sizeof(Dtype), 1, file);
        }
      }
    }
    fclose(file);

    hipMemcpy(fieldCI + n * nx * ny * nz, ch[n].fieldCI, sizeof(Dtype) * nx * ny * nz, hipMemcpyHostToDevice);
  }
#endif
}
#if 1
void
write_field2B (int irun, int mode)
{
	int i, j, k, n;
	char filename[1024];
	FILE *file;
	if(prank == myrank){
		printf("write_chk %d\n", irun);
	}

	switch (mode) {
	case 0:
		for(n = 0; n < nac; n++)
		{
			hipMemcpy (ac[n].field2b, field2BE + n * offset, sizeof (Stype) * offset, hipMemcpyDeviceToHost);

			sprintf (filename, "%s%02d/eta%d_%06d_%02d%02d%02d.dat", data_dir, cart_id[2], n, irun, cart_id[0], cart_id[1], cart_id[2]);
			file = fopen (filename, "wb");
			for (k = iz1; k < iz4; k++) {
				for (j = iy1; j < iy4; j++) {
					for (i = ix1; i < ix4; i++) {
						fwrite (&ac[n].field2b[k * nx * ny + j * nx + i], sizeof (Stype), 1, file);
					}
				}
			}
			fclose (file);
		}
		break;
	case 1:
    #if 0
		for(n = 0; n < nch; n++)
		{
			hipMemcpy (ch[n].field2b, field2BCI + n * offset, sizeof (Stype) * offset, hipMemcpyDeviceToHost);

			sprintf (filename, "%s%02d/c%d_%06d_%02d%02d%02d.dat", data_dir, cart_id[2], n, irun, cart_id[0], cart_id[1], cart_id[2]);
			file = fopen (filename, "wb");
			for (k = iz1; k < iz4; k++) {
				for (j = iy1; j < iy4; j++) {
					for (i = ix1; i < ix4; i++) {
						fwrite (&ch[n].field2b[k * nx * ny + j * nx + i], sizeof (Stype), 1, file);
					}
				}
			}
			fclose (file);
		}
    #endif
		break;
 	}
}
#endif
