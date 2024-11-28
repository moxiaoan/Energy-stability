#include   <stdio.h>
#include   "mpi.h"
#include   "ScLETD.h"
#include   <string.h>
#include   <errno.h>
#include   <stdlib.h>


int   main(int   argc,   char   **argv)
{
      int   n;
      double   tt;
      char   filename[1024];
      FILE   *file;
      start_mpi(argc,   argv);
      sprintf(filename,   "input%s.ini",   argv[1]);
      file   =   fopen(filename,   "r");
      fscanf(file,   "%d,%d,%d",   &nx,   &ny,   &nz);
      fscanf(file,   "%d,%d,%d",   &procs[0],   &procs[1],   &procs[2]);
      fscanf(file,   "%d",   &nghost);
      fscanf(file,   "%s",   work_dir);
      fscanf(file,   "%s",   data_dir);
      fscanf(file,   "%lf,%lf",   &dt,   &t_total);
      fscanf(file,   "%lf",   &epn2);
      fscanf(file,   "%d",   &periodic);
      fscanf(file,   "%d,   %d",   &nac,   &nch);
      fscanf(file,   "%d,   %d,   %s",   &restart,   &restart_iter,   restart_dir);
      fscanf(file,   "%lf,%lf,%lf",   &xmin,   &ymin,   &zmin);
      fscanf(file,   "%lf,%lf,%lf",   &xmax,   &ymax,   &zmax);
      fscanf(file,   "%d",   &ETD2);
      fscanf(file,   "%d",   &ELASTIC);
      fscanf(file,   "%d",   &ANISOTROPIC);
      fscanf(file,   "%lf,%lf,%lf",   &kkx,   &kky,   &kkz);
      fscanf(file,   "%d",   &VariableMobility);
      fclose(file);
      if   (myrank   ==   prank)
      {
            printf("running   on\t%d   processors\n",   nprocs);
            printf("nx,ny,nz\t%d\t%d\t%d\n",   nx,   ny,   nz);
            printf("px,py,pz\t%d\t%d\t%d\n",   procs[0],   procs[1],   procs[2]);
            printf("nghost\t\t%d\n",   nghost);
            printf("ini   file\t%s\n",   filename);
            printf("output   dir\t%s\n",   work_dir);
            printf("data   dir\t%s\n",   data_dir);
            printf("epn2\t\t%lf\n",   epn2);
            printf("dt\t\t%lf\n",   dt);
            printf("nac,nch\t\t%d,%d\n",   nac,   nch);
            printf("restart,   iter,   dir\t%d,%d,%s\n",   restart,   restart_iter,   restart_dir);
            printf("ETD2\t\t%d\n",   ETD2);
            printf("ELASTIC\t\t%d\n",   ELASTIC);
            printf("ANISOTROPIC\t\t%d\n",   ANISOTROPIC);
            printf("VariableMobility\t\t%d\n",   VariableMobility);
      }
      if   (nac   ==   0   &&   nch   ==   0)
      {
            printf("nac   =   0   and   nch   =   0\n");
            exit(1);
      }
      if   (ELASTIC   !=   NEITHER_FUNCTION)
      {
            elastic_input();
      }
      alloc_vars();
      init_vars();
      if   (ANISOTROPIC   ==   1)
      {
            anisotropic_input();
      }
      define_mpi_type();
      sw_cart_creat();
      init_para();
      read_matrices();
      init_field();
      if   (ELASTIC   !=   NEITHER_FUNCTION)
      {
            elastic_init();
      }
      iter   =   0;
      if   (restart   ==   1)   {
            iter   =   restart_iter;
      }
      tt   =   iter   *   dt;
      ioutput   =   iter   /   100;
      while   (tt   <   t_total)
      {
            check_soln_new(tt);
            iter++;
            tt   +=   dt;
            transfer();
            if   (ELASTIC   !=   NEITHER_FUNCTION)
            {
                  elastic_calculate();
            }
            calc_mu();

            for   (n   =   0;   n   <   nac;   n++)
            {
                  ac_calc_FU(n,   ac[n].fieldE1);
            }
            for   (n   =   0;   n   <   nch;   n++)
            {
                  ch_calc_FU(n,   ch[n].fieldCI1);
            }

            for   (n   =   0;   n   <   nac;   n++)
            {
                  stage   =   0;
                  PUX(MPXI,   ac[n].fieldE,   ac[n].fieldEt,   ac[n].fieldEp);
                  PUY(MPYI,   ac[n].fieldEt,   ac[n].fieldE,   ac[n].fieldEp);
                  PUZ(MPZI,   ac[n].fieldE,   ac[n].fieldEt,   ac[n].fieldEp,   ac[n].fieldEt);

                  stage   =   1;
                  PUX(MPXI,   ac[n].fieldE1,   ac[n].fieldE,   ac[n].fieldEp);
                  PUY(MPYI,   ac[n].fieldE,   ac[n].fieldE1p,   ac[n].fieldEp);
                  PUZ(MPZI,   ac[n].fieldE1p,   ac[n].fieldE,   ac[n].fieldEp,   ac[n].fieldEt);

                  ac_updateU_new(n,   ac[n].fieldEt,   ac[n].fieldEp);
                  zxy_xyz(ac[n].fieldEt,   ac[n].fieldE);

                  stage   =   2;
                  PUX(MPX,   ac[n].fieldE,   ac[n].fieldE1p,   ac[n].fieldEp);
                  PUY(MPY,   ac[n].fieldE1p,   ac[n].fieldEt,   ac[n].fieldEp);
                  PUZ(MPZ,   ac[n].fieldEt,   ac[n].fieldE,   ac[n].fieldEp,   ac[n].fieldEt);
            }

            for   (n   =   0;   n   <   nch;   n++)
            {
                  stage   =   0;
                  PUX(MPXI,   ch[n].fieldCI,   ch[n].fieldCIt,   ch[n].fieldCIp);
                  PUY(MPYI,   ch[n].fieldCIt,   ch[n].fieldCI,   ch[n].fieldCIp);
                  PUZ(MPZI,   ch[n].fieldCI,   ch[n].fieldCIt,   ch[n].fieldCIp,   ch[n].fieldCIt);
                  stage   =   1;
                  PUX(MPXI,   ch[n].fieldCI1,   ch[n].fieldCI,   ch[n].fieldCIp);
                  PUY(MPYI,   ch[n].fieldCI,   ch[n].fieldCI1p,   ch[n].fieldCIp);
                  PUZ(MPZI,   ch[n].fieldCI1p,   ch[n].fieldCI,   ch[n].fieldCIp,   ch[n].fieldCIt);
                  ch_updateU_new(n,   ch[n].fieldCIt,   ch[n].fieldCIp);
                  zxy_xyz(ch[n].fieldCIt,   ch[n].fieldCI);
                  stage   =   2;
                  PUX(MPX,   ch[n].fieldCI,   ch[n].fieldCI1p,   ch[n].fieldCIp);
                  PUY(MPY,   ch[n].fieldCI1p,   ch[n].fieldCIt,   ch[n].fieldCIp);
                  PUZ(MPZ,   ch[n].fieldCIt,   ch[n].fieldCI,   ch[n].fieldCIp,   ch[n].fieldCIt);
            }
            if   (ETD2   ==   1)
            {
                  transfer();
                  if   (ELASTIC   !=   NEITHER_FUNCTION)
                  {
                        elastic_calculate();
                  }
                  calc_mu();
                  for   (n   =   0;   n   <   nac;   n++)
                  {
                        ac_calc_FU(n,   ac[n].fieldE2);
                  }
                  for   (n   =   0;   n   <   nch;   n++)
                  {
                        ch_calc_FU(n,   ch[n].fieldCI2);
                  }
                  for   (n   =   0;   n   <   nac;   n++)
                  {
                        prepare_U1_new(ac[n].fieldE1,   ac[n].fieldE2);
                        stage   =   0;
                        PUX(MPXI,   ac[n].fieldE2,   ac[n].fieldEt,   ac[n].fieldEp);
                        PUY(MPYI,   ac[n].fieldEt,   ac[n].fieldE1p,   ac[n].fieldEp);
                        PUZ(MPZI,   ac[n].fieldE1p,   ac[n].fieldE2,   ac[n].fieldEp,   ac[n].fieldEt);
                        prepare_U2_new(ac[n].phiE2,   ac[n].fieldE1,   ac[n].fieldE2);
                        zxy_xyz(ac[n].fieldE2,   ac[n].fieldE1);
                        stage   =   2;
                        PUX(MPX,   ac[n].fieldE1,   ac[n].fieldE1p,   ac[n].fieldEp);
                        PUY(MPY,   ac[n].fieldE1p,   ac[n].fieldEt,   ac[n].fieldEp);
                        PUZ(MPZ,   ac[n].fieldEt,   ac[n].fieldE1,   ac[n].fieldEp,   ac[n].fieldEt);
                        correct_U_new(ac[n].fieldE,   ac[n].fieldE1);
                  }
                  for   (n   =   0;   n   <   nch;   n++)
                  {
                        prepare_U1_new(ch[n].fieldCI1,   ch[n].fieldCI2);
                        stage   =   0;
                        PUX(MPXI,   ch[n].fieldCI2,   ch[n].fieldCIt,   ch[n].fieldCIp);
                        PUY(MPYI,   ch[n].fieldCIt,   ch[n].fieldCI1p,   ch[n].fieldCIp);
                        PUZ(MPZI,   ch[n].fieldCI1p,   ch[n].fieldCI2,   ch[n].fieldCIp,   ch[n].fieldCIt);
                        prepare_U2_new(ch[n].phiCI2,   ch[n].fieldCI1,   ch[n].fieldCI2);
                        zxy_xyz(ch[n].fieldCI2,   ch[n].fieldCI1);
                        stage   =   2;
                        PUX(MPX,   ch[n].fieldCI1,   ch[n].fieldCI1p,   ch[n].fieldCIp);
                        PUY(MPY,   ch[n].fieldCI1p,   ch[n].fieldCIt,   ch[n].fieldCIp);
                        PUZ(MPZ,   ch[n].fieldCIt,   ch[n].fieldCI1,   ch[n].fieldCIp,   ch[n].fieldCIt);
                        correct_U_new(ch[n].fieldCI,   ch[n].fieldCI1);
                  }
            }
            MPI_Barrier(MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      if   (ELASTIC   !=   NEITHER_FUNCTION)
      {
            elastic_finish();
      }
      dealloc_vars();
      close_mpi();

      return   0;
}
