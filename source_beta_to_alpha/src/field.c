#include "ScLETD.h"
#include<sys/types.h>
#include<unistd.h>
#include<limits.h>
#include<sys/stat.h>
#include "anisotropic_hip.h"

void
outputfield ()
{
  int i, j, k, n;
  char filename[1024];
  FILE *file;

for (n = 0; n < nac; n++) 
{
  //sprintf (filename, "%seta%d_%06d_%02d%02d%02d.dat", data_dir, n, ioutput, cart_id[0], cart_id[1], cart_id[2]);
  sprintf (filename, "%s%02d/eta%d_%06d_%02d%02d%02d.dat", data_dir, cart_id[2], n, ioutput, cart_id[0], cart_id[1], cart_id[2]);
  //file = fopen (filename, "wb");
  file = fopen (filename, "w");
  //if (file == NULL) {
   // printf ("fopen error %s!\n", strerror(errno));
   // exit (1);
  //}
  for (k = nghost; k < nz - nghost; k++) {
    for (j = nghost; j < ny - nghost; j++) {
      for (i = nghost; i < nx - nghost; i++) {
        //fwrite (&ac[n].fieldE[k * nx * ny + j * nx + i], sizeof (double), 1, file);
        fprintf (file, "%+1.15lf\n", ac[n].fieldE[k * nx * ny + j * nx + i]);
      }
    }
  }
  fclose (file);
}
for (n = 0; n < nch; n++) 
{
  //sprintf (filename, "%sc%d_%06d_%02d%02d%02d.dat", data_dir, n, ioutput, cart_id[0], cart_id[1], cart_id[2]);
  sprintf (filename, "%s%02d/c%d_%06d_%02d%02d%02d.dat", data_dir, cart_id[2], n, ioutput, cart_id[0], cart_id[1], cart_id[2]);
 // file = fopen (filename, "wb");
  file = fopen (filename, "w");

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        //fwrite (&ch[n].fieldCI[k * nx * ny + j * nx + i], sizeof (double), 1, file);
        fprintf (file, "%+1.15lf\n", ch[n].fieldCI[k * nx * ny + j * nx + i]);
      }
    }
  }
  fclose (file);
}

  ioutput+=1;
}

void
write_section (char *fname, double *field)
{
  FILE *fp;
  int i, j, k;
  if (cart_id[2] == (procs[2]/2)) {
    fp = fopen (fname, "w");
    if (fp == NULL) {
      printf ("fopen error %s!\n", strerror(errno));
      exit (1);
    }
    if (procs[2] == 1) {
      k = nz / 2;
    } else {
      k = 2;
    }
    for (j = nghost; j < ny-nghost; j++) {
      for (i = nghost; i < nx-nghost; i++) {
        fprintf (fp, "%+1.15lf ", field[k * nx * ny + j * nx + i]);
      }
      fprintf (fp, "\n");
    }
    fclose (fp);
  }
}

void
read_ori (void)
{
  FILE *fp;
  int i, j, k;
  char filename[1024];

  //sprintf (filename, "/work1/zhangjian/alpha_beta_test/1024_node_run/TEST/eta12_c2_dgemm_elastic_GX_noelastic_AI/data/%02d/0%02d_0%02d_0%02d", cart_id[2], cart_id[0]+1, cart_id[1]+1, cart_id[2]+1);
  sprintf (filename, "/work1/zhangjian/alpha_beta_test/1024_node_run/tmp_444/0%02d_0%02d_0%02d", cart_id[0]+1, cart_id[1]+1, cart_id[2]+1);
  //sprintf (filename, "./data/%02d/ori_7_%02d%02d%02d.dat", cart_id[2], cart_id[0], cart_id[1], cart_id[2]);
  //fp = fopen (filename, "r");
  fp = fopen (filename, "rb");
  if (fp == NULL) {
    printf ("fopen error %s!\n", strerror(errno));
    exit (1);
  }
// for (i = 0; i < nx; i++) {
//   for (j = 0; j < ny; j++) {
//     for (k = 0; k < nz; k++) {
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
//        fscanf (fp, "%d", &ori[k * nx * ny + j * nx + i]);
	fread (&ori[k * nx * ny + j * nx + i], sizeof (int), 1, fp);
      }
    }
  }
  fclose (fp);
  hipMemcpy (Ori, ori, nx * ny * nz * sizeof (int), hipMemcpyHostToDevice);
#if 1
  sprintf (filename, "/work1/zhangjian/alpha_beta_test/1024_node_run/tmp_444/f_7_%02d%02d%02d.dat", cart_id[0], cart_id[1], cart_id[2]);
  //sprintf (filename, "./data/%02d/f_7_%02d%02d%02d.dat", cart_id[2], cart_id[0], cart_id[1], cart_id[2]);
  //fp = fopen (filename, "r");
  fp = fopen (filename, "rb");
  if (fp == NULL) {
    printf ("fopen error %s!\n", strerror(errno));
    exit (1);
  }
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
//        fscanf (fp, "%d", &f[k * nx * ny + j * nx + i]);
	fread (&f[k * nx * ny + j * nx + i], sizeof (int), 1, fp);
      }
    }
  }
  fclose (fp);
  hipMemcpy (F, f, nx * ny * nz * sizeof (int), hipMemcpyHostToDevice);
#endif
}

void
read_field (char *fname, double *field)
{
  FILE *fp;
  int i, j, k;
  fp = fopen (fname, "r");
  if (fp == NULL) {
    printf ("fopen error %s!\n", strerror(errno));
    exit (1);
  }
  for (k = nghost; k < nz - nghost; k++) {
    for (j = nghost; j < ny - nghost; j++) {
      for (i = nghost; i < nx - nghost; i++) {
        fscanf (fp, "%lf", &field[k * nx * ny + j * nx + i]);
      }
    }
  }
  fclose (fp);
}

void
write_field (char* fname, double* field)
{
  FILE *fp;
  int i, j, k;
  fp = fopen (fname, "w");
  if (fp == NULL) {
    printf ("fopen error %s!\n", strerror(errno));
    exit (1);
  }
  for (k = nghost; k < nz-nghost; k++) {
    for (j = nghost; j < ny-nghost; j++) {
      for (i = nghost; i < nx-nghost; i++) {
        fprintf (fp, "%+1.15lf\n", field[k * nx * ny + j * nx + i]);
      }
    }
  }
  fclose (fp);
}


void print_info(void)
{
     int buf[1];
     MPI_Status status1;
     buf[0] = 1;
     if ( myrank>0) {
         MPI_Recv(buf, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, &status1);
     }
      
     printf("step%d, myrank = %d, (%d, %d, %d), on %s\n",iter,myrank,cart_id[0],cart_id[1],cart_id[2],processor_name);
 
     if ( myrank<nprocs-1) {
         MPI_Send(buf, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD);
     }
     MPI_Barrier(MPI_COMM_WORLD);
}

void
//check_cahn_hilliard_min (double *fieldci0_old, double *fieldci1_old, int flag)
check_cahn_hilliard_min (int &flag)
{
	int i, j, k;
	double tmp, mintmp, mmineta0, mmineta1, mminc0, mminc1;
#if 0
	mintmp = 1.0e30;
	for (k = iz2; k < iz3; k++) {
		for (j = iy2; j < iy3; j++) {
			for (i = ix2; i < ix3; i++) {
				tmp = ac[0].fieldE[k * nx * ny + j * nx + i];
				if (tmp < mintmp) {
					mintmp = tmp;
				}
			}
		}
	}
	MPI_Reduce (&mintmp, &mmineta0, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);
	mintmp = 1.0e30;
	for (k = iz2; k < iz3; k++) {
		for (j = iy2; j < iy3; j++) {
			for (i = ix2; i < ix3; i++) {
				tmp = ac[1].fieldE[k * nx * ny + j * nx + i];
				if (tmp < mintmp) {
					mintmp = tmp;
				}
			}
		}
	}
	MPI_Reduce (&mintmp, &mmineta1, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);
#endif
	mintmp = 1.0e30;
	for (k = iz2; k < iz3; k++) {
		for (j = iy2; j < iy3; j++) {
			for (i = ix2; i < ix3; i++) {
				tmp = ch[0].fieldCI[k * nx * ny + j * nx + i];
				if (tmp < mintmp) {
					mintmp = tmp;
				}
			}
		}
	}
	MPI_Reduce (&mintmp, &mminc0, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);
	mintmp = 1.0e30;
	for (k = iz2; k < iz3; k++) {
		for (j = iy2; j < iy3; j++) {
			for (i = ix2; i < ix3; i++) {
				tmp = ch[1].fieldCI[k * nx * ny + j * nx + i];
				if (tmp < mintmp) {
					mintmp = tmp;
				}
			}
		}
	}
	MPI_Reduce (&mintmp, &mminc1, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);
	if (myrank == prank) {
		//if ((mmin0 <= 0.049) || (mmin1 <= 0.009))
		//if (((mminc0 <= 0.04) || (mminc1 <= 0.007)) || ((mmineta0 <= -0.02) || (mmineta0 <= -0.02)))
		//if ((mminc0 <= 0.04) || (mminc1 <= 0.007))
		//if ((mminc0 <= 0.049) || (mminc1 <= 0.009))
		if ((mminc0 <= 0.039) || (mminc1 <= 0.007))
		//if ((mminc0 <= 0.025) || (mminc1 <= 0.003))
		{
			flag = 1;
		}
		//if ((mminc0 >= 0.05) && (mminc1 >= 0.0095) && (mmineta0 >= 0.0) && (mmineta0 >= 0.0))
		if ((mminc0 >= 0.065) && (mminc1 >= 0.0095))
		//if ((mminc0 >= 0.045) && (mminc1 >= 0.0075))
		{
			flag = 2;
		}
	}
	MPI_Bcast (&flag, 1, MPI_INT, prank, MPI_COMM_WORLD);
}

void
check_soln_new (double time)
{
  if (myrank == prank) {
    printf ("--------------iter %d----count %d----------\n", iter, count_gyq);
  }
  int n;
  int i, j, k;
  double tmp, maxtmp, mintmp, voltmp, mmax, mmin, vol;
  FILE *fp;
  char fname[1024];
  if (iter == 0) {
    for (n = 0; n < nac; n++) {
      if (myrank == prank) {
        sprintf (fname, "%sphi_ac_%d.csv", work_dir, n);
        fp = fopen (fname, "w");
        if (fp == NULL) {
          printf ("fopen error %s!\n", strerror(errno));
          exit (1);
        }
        fprintf (fp, "iter,maxphiE,minphiE,volE\n");
        fclose (fp);
      }
    }
    for (n = 0; n < nch; n++) {
      if (myrank == prank) {
        sprintf (fname, "%sphi_ch_%d.csv", work_dir, n);
        fp = fopen (fname, "w");
        if (fp == NULL) {
          printf ("fopen error %s!\n", strerror(errno));
          exit (1);
        }
        fprintf (fp, "iter,maxphiCI,minphiCI,volCI\n");
        fclose (fp);
      }
    }
  }
  for (n = 0; n < nac; n++) {
    if (myrank == prank) {
      sprintf (fname, "%sphi_ac_%d.csv", work_dir, n);
      errno = 0;
      fp = fopen (fname, "a");
      if (fp == NULL) {
        printf ("fopen error %s!\n", strerror(errno));
        exit (1);
      }
    }
    maxtmp = -1.0e30;
    mintmp = 1.0e30;
    voltmp = 0.0;
    for (k = iz2; k < iz3; k++) {
      for (j = iy2; j < iy3; j++) {
        for (i = ix2; i < ix3; i++) {
           tmp = ac[n].fieldE[k * nx * ny + j * nx + i];
          if (tmp > maxtmp) {
            maxtmp = tmp;
          }
          if (tmp < mintmp) {
            mintmp = tmp;
          }
          voltmp += tmp + 1.0;
        }
      }
    }
    voltmp *= hx * hy * hz;
    #if 0
    double c0_bcc = 0.102;
    double c1_bcc = 0.036;
    double c2_bcc = 1.0 - c0_bcc - c1_bcc;
    double c0_hcp = 0.1057;
    double c1_hcp = 0.0223;
    double c2_hcp = 1.0 - c0_hcp - c1_hcp;
  double GBCC = c2_bcc * GBCCTI + c0_bcc * GBCCAL + c1_bcc * GHSERV + R_a * T_a * (c2_bcc * log(c2_bcc) + c0_bcc * log(c0_bcc) + c1_bcc * log(c1_bcc)) + c0_bcc * c1_bcc * (BL12_0 + BL12_1 * (c0_bcc - c1_bcc)) + c2_bcc * c1_bcc * (BL32_0 + BL32_1 * (c2_bcc - c1_bcc) + BL32_2 * (c2_bcc - c1_bcc) * (c2_bcc - c1_bcc)) + c0_bcc * c2_bcc * (BL13_0 + BL13_1 * (c0_bcc - c2_bcc) + BL13_2 * (c0_bcc - c2_bcc) * (c0_bcc - c2_bcc)) + c0_bcc * c1_bcc * c2_bcc * BL132_0; 

  double GHCP = c2_hcp * GHSERTI + c0_hcp * GHCPAL + c1_hcp * GHCPV + R_a * T_a * (c2_hcp * log(c2_hcp) + c0_hcp * log(c0_hcp) + c1_hcp * log(c1_hcp)) + c0_hcp * c1_hcp * (HL12_0 + HL12_1 * (c0_hcp - c1_hcp)) + c2_hcp * c1_hcp * HL32_0 + c0_hcp * c2_hcp * (HL13_0 + HL13_1 * (c0_hcp - c2_hcp) + HL13_2 * (c0_hcp - c2_hcp) * (c0_hcp - c2_hcp)) + c0_hcp * c1_hcp * c2_hcp * (c0_hcp * HL132_0 + c2_hcp * HL132_1 + c1_hcp * HL132_2);
  printf("GBCC = %lf, GHCP = %lf\n", GBCC, GHCP);
  #endif
    MPI_Reduce (&maxtmp, &mmax, 1, MPI_DOUBLE, MPI_MAX, prank, MPI_COMM_WORLD);
    MPI_Reduce (&mintmp, &mmin, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);
    MPI_Reduce (&voltmp, &vol, 1, MPI_DOUBLE, MPI_SUM, prank, MPI_COMM_WORLD);
//    printf("step%dfieldE, myrank = %d, (%d, %d, %d), on %s\n",iter,myrank,cart_id[0],cart_id[1],cart_id[2],processor_name);
    if (myrank == prank) {
      printf ("ac[%d].fieldE\n", n);
      printf ("max phi\t\t%+1.15lf\n", mmax);
      printf ("min phi\t\t%+1.15lf\n", mmin);
      printf ("vol\t\t%+1.15lf\n", vol);
      fprintf (fp, "%d,%+1.15lf,%+1.15lf,%+1.15lf\n", iter, mmax, mmin, vol);
      fclose (fp); 
    }
  }
#if 0
  for (n = 0; n < nch; n++) {
    if (myrank == prank) {
      sprintf (fname, "%sphi_ch_%d.csv", work_dir, n);
      fp = fopen (fname, "a");
      if (fp == NULL) {
        printf ("fopen error %s!\n", strerror(errno));
        exit (1);
      }
    }  
    maxtmp = -1.0e30;
    mintmp = 1.0e30;
    voltmp = 0.0;

    for (k = iz2; k < iz3; k++) {
      for (j = iy2; j < iy3; j++) {
        for (i = ix2; i < ix3; i++) {
          tmp = ch[n].fieldCI[k * nx * ny + j * nx + i];
          if (tmp > maxtmp) {
            maxtmp = tmp;
          }
          if (tmp < mintmp) {
            mintmp = tmp;
          }
          voltmp += tmp + 1;
        }
      }
    }
    voltmp *= hx * hy * hz;
    MPI_Reduce (&maxtmp, &mmax, 1, MPI_DOUBLE, MPI_MAX, prank, MPI_COMM_WORLD);
    MPI_Reduce (&mintmp, &mmin, 1, MPI_DOUBLE, MPI_MIN, prank, MPI_COMM_WORLD);
    MPI_Reduce (&voltmp, &vol, 1, MPI_DOUBLE, MPI_SUM, prank, MPI_COMM_WORLD);

    if (myrank == prank) {
      printf ("ch[%d].fieldCI\n", n);
      printf ("max phi\t\t%+1.15lf\n", mmax);
      printf ("min phi\t\t%+1.15lf\n", mmin);
      printf ("vol\t\t%+1.15lf\n", vol);
      fprintf (fp, "%d,%+1.15lf,%+1.15lf,%+1.15lf\n", iter, mmax, mmin, vol);
      fclose (fp); 
    }
  }
#endif
//  if (iter != 0 && iter % 5000 == 0)
//    outputfield ();
}
/*
void
init_field_origin() {
  int i, j, k;
  int gx, gx1;
  int pos1, pos2;
  float HARVEST1, HARVEST2;
  
  for (k = iz1; k < iz4; k++)
  {
     for (j = iy1; j < iy4; j++)
     {
       HARVEST2 = rand() / (double)(RAND_MAX);
       pos2 = (procs[0] * nx - 20) + 3.0 * (HARVEST1 - 0.5);
       HARVEST1 = rand() / (double)(RAND_MAX);
       pos1 = 19 + 3.0 * (HARVEST2 - 0.5);
       for (i = ix1; i < ix4; i++)
       {
       	   gx = cart_id[0] * nx + i;
       	   //gx1 = (cart_id[0] + 1) * nx - 1 + i;
	   if (gx <= pos1)
	   {
	   	ac[0].fieldE[k * ny * nx + j * nx + i] = 1.0;	
	   	ac[1].fieldE[k * ny * nx + j * nx + i] = 0.0;	
	        ch[0].fieldCI[k * ny * nx + j * nx + i] = 0.1057;
	        ch[1].fieldCI[k * ny * nx + j * nx + i] = 0.0223;
	   }
	   else if (gx >= pos2)
	   {
	   	ac[0].fieldE[k * ny * nx + j * nx + i] = 0.0;	
	   	ac[1].fieldE[k * ny * nx + j * nx + i] = 1.0;	
	        ch[0].fieldCI[k * ny * nx + j * nx + i] = 0.1057;
	        ch[1].fieldCI[k * ny * nx + j * nx + i] = 0.0223;
	   }
	   else
	   {
	   	ac[0].fieldE[k * ny * nx + j * nx + i] = 0.0;	
	   	ac[1].fieldE[k * ny * nx + j * nx + i] = 0.0;	
	     	ch[0].fieldCI[k * ny * nx + j * nx + i] = 0.102;
	     	ch[1].fieldCI[k * ny * nx + j * nx + i] = 0.036;
	   }
       }
     }
  }
  for (int m = 0; m < nac; m++)
  {
     hipMemcpy (fieldE + m * offset, ac[m].fieldE, nx * ny * nz * sizeof (double), hipMemcpyHostToDevice);
  }
  for (int m = 0; m < nch; m++)
  {
     hipMemcpy (fieldCI + m * offset, ch[m].fieldCI, nx * ny * nz * sizeof (double), hipMemcpyHostToDevice);
  }
}*/

void
init_field_origin() {
  int i, j, k;
  double rad = 10.0;
  
  for (k = iz1; k < iz4; k++)
  {
     for (j = iy1; j < iy4; j++)
     {
       for (i = ix1; i < ix4; i++)
       {
	   if ((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<rad*rad)
	   {
	     ac[0].fieldE[k * ny * nx + j * nx + i] = 1.0;
	     //ac[1].fieldE[k * ny * nx + j * nx + i] = 0.0;
	   }
	   else
	   {
	     ac[0].fieldE[k * ny * nx + j * nx + i] = 0.0;
	   }
       }
     }
  }
}
/*
void initial_ac (double *fielde1, double *fielde2, int n)
{
  double thr;
  double rad = 3.0;
  int i, j, k;
  int x, y, z;
  thr = rand() / (double)(RAND_MAX);
  if (n == 0)
  {
  if (thr < 0.9)
  {
	x = (int)((rand() / (double)(RAND_MAX)) * nx);
	y = (int)((rand() / (double)(RAND_MAX)) * ny);
	z = (int)((rand() / (double)(RAND_MAX)) * nz);

	for (i = ix1 + rad; i < ix4 - rad; i++)
	{
		for (j = iy1 + rad; j < iy4 - rad; j++)
		{
			for (k = iz1 + rad; k < iz4 - rad; k++)
			{
				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad)
				{
					fielde1[k * nx * ny + j * nx + i] = 1.0;
				}				
//			x = abs(x-i);
//			y = abs(y-j);
//			z = abs(z-k);
//			if ((x<rad) && (y<rad) && (z<rad)){
//				fielde1[k * nx * ny + j * nx + i] = 1.0;
//			}
			}
		}
        
	}

  }
  }
  else
  {
  if (thr < 0.9)
  {
	x = (int)((rand() / (double)(RAND_MAX)) * nx);
	y = (int)((rand() / (double)(RAND_MAX)) * ny);
	z = (int)((rand() / (double)(RAND_MAX)) * nz);

	for (i = ix1 + rad; i < ix4 - rad; i++)
	{
		for (j = iy1 + rad; j < iy4 - rad; j++)
		{
			for (k = iz1 + rad; k < iz4 - rad; k++)
			{
				int l = k * nx * ny + j * nx + i;
				if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && 					fielde1[l] == 0.0)
				{
					fielde2[k * nx * ny + j * nx + i] = 1.0;
				}				
//			x = abs(x-i);
//			y = abs(y-j);
//			z = abs(z-k);
//			if ((x<rad) && (y<rad) && (z<rad) && (fielde1[l] == 0.0)){
//				fielde2[k * nx * ny + j * nx + i] = 1.0;
//			}
			}
		}
        
	}

  }

  }
}*/

void initial_zero (double *fielde1, double *fielde2, double *fielde3, double *fielde4, double *fielde5, double *fielde6)
{
  int i, j, k;
  for (i = ix1; i < ix4; i++)
  {
  	for (j = iy1; j < iy4; j++)
  	{
         	for (k = iz1; k < iz4; k++)
           	{
           		fielde1[k * nx * ny + j * nx + i] = 0.0;
           		fielde2[k * nx * ny + j * nx + i] = 0.0;
           		fielde3[k * nx * ny + j * nx + i] = 0.0;
           		fielde4[k * nx * ny + j * nx + i] = 0.0;
           		fielde5[k * nx * ny + j * nx + i] = 0.0;
           		fielde6[k * nx * ny + j * nx + i] = 0.0;
		    }
	}
  }
}
/*
void initial_ac (double *fielde1, double *fielde2, double *fielde3, double *fielde4, double *fielde5, double *fielde6, double *fielde7, double *fielde8, double *fielde9, double *fielde10, double *fielde11, double *fielde12, int n)
{
  double thr;
  double rad = 5.0;
  int i, j, k;
  int x, y, z;
  thr = rand() / (double)(RAND_MAX);
  if (n == 0)
  {
        if (thr < 0.65)
        {
              if (cart_id[0] == 0 && cart_id[1] == 0)
              {
               	      x = (int)((rand() / (double)(RAND_MAX)) * nx);
              	      y = (int)((rand() / (double)(RAND_MAX)) * ny);
              	      z = (int)((rand() / (double)(RAND_MAX)) * nz);
              
                      x = x + ((rad - x) * (x <= rad));
                      y = y + ((rad - x) * (y <= rad));
                      z = z + ((rad - x) * (z <= rad));
              
                      x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
                      y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
                      z = z - ((z - (nz - rad)) * (z >= (nz - rad)));
              
              
              	      for (i = ix1 + rad; i < ix4 - rad; i++)
              	      {
              	       	    for (j = iy1 + rad; j < iy4 - rad; j++)
              		    {
              			  for (k = iz1 + rad; k < iz4 - rad; k++)
              			  {
              				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad)
              				{
              					fielde1[k * nx * ny + j * nx + i] = 1.0;
              				}				
              			  }
              		    }
              	      }
              }
        }
  }
  else if (n == 1)
  {
	if (thr < 0.65)
	{
		if (cart_id[0] == 0 && cart_id[1] == 0)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde1[l] == 0.0)
						{
							fielde2[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 2)
  {
	if (thr < 0.65)
	{
		if (cart_id[0] == 0 && cart_id[1] == 0)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde1[l] == 0.0 && fielde2[l] == 0.0)
						{
							fielde3[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 3)
  {
	if (thr < 0.75)
	{
		if (cart_id[0] == 1 && cart_id[1] == 0)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad))
						{
							fielde4[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 4)
  {
	if (thr < 0.75)
	{
		if (cart_id[0] == 1 && cart_id[1] == 0)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde4[l] == 0.0)
						{
							fielde5[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 5)
  {
	if (thr < 0.55)
	{
		if (cart_id[0] == 1 && cart_id[1] == 0)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde4[l] == 0.0 && fielde5[l] == 0.0)
						{
							fielde6[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 6)
  {
	if (thr < 0.65)
	{
		if (cart_id[0] == 0 && cart_id[1] == 1)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad))
						{
							fielde7[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 7)
  {
	if (thr < 0.55)
	{
		if (cart_id[0] == 0 && cart_id[1] == 1)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde7[l] == 0.0)
						{
							fielde8[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 8)
  {
	if (thr < 0.65)
	{
		if (cart_id[0] == 0 && cart_id[1] == 1)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde7[l] == 0.0 && fielde8[l] == 0.0)
						{
							fielde9[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 9)
  {
	if (thr < 0.75)
	{
		if (cart_id[0] == 1 && cart_id[1] == 1)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad))
						{
							fielde10[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else if (n == 10)
  {
	if (thr < 0.55)
	{
		if (cart_id[0] == 1 && cart_id[1] == 1)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde10[l] == 0.0)
						{
							fielde11[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
  	}
  }
  else
  {
	if (thr < 0.55)
	{
		if (cart_id[0] == 1 && cart_id[1] == 1)
		{
			x = (int)((rand() / (double)(RAND_MAX)) * nx);
			y = (int)((rand() / (double)(RAND_MAX)) * ny);
			z = (int)((rand() / (double)(RAND_MAX)) * nz);

        		x = x + ((rad - x) * (x <= rad));
        		y = y + ((rad - x) * (y <= rad));
        		z = z + ((rad - x) * (z <= rad));

        		x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        		y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        		z = z - ((z - (nz - rad)) * (z >= (nz - rad)));

			for (i = ix1 + rad; i < ix4 - rad; i++)
			{
				for (j = iy1 + rad; j < iy4 - rad; j++)
				{
					for (k = iz1 + rad; k < iz4 - rad; k++)
					{
						int l = k * nx * ny + j * nx + i;
						if ((sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad) && fielde10[l] == 0.0 && fielde11[l] == 0.0)
						{
							fielde12[k * nx * ny + j * nx + i] = 1.0;
						}				
					}
				}
			}
		}
	}
  }
}*/
/*
void initial_ac (int *orie, double *fielde1, double *fielde2, double *fielde3, double *fielde4, double *fielde5, double *fielde6, double *fielde7, double *fielde8, double *fielde9, double *fielde10, double *fielde11, double *fielde12)
{
  double thr;
  double rad = 5.0;
  int i, j, k;
  int x, y, z, l;
  thr = rand() / (double)(RAND_MAX);
  int o; 

  if (thr < 0.9)
  {
         x = (int)((rand() / (double)(RAND_MAX)) * nx);
         y = (int)((rand() / (double)(RAND_MAX)) * ny);
         z = (int)((rand() / (double)(RAND_MAX)) * nz);
        
        //x = x + ((rad - x) * (x <= rad));
        //y = y + ((rad - x) * (y <= rad));
        //z = z + ((rad - x) * (z <= rad));
        
        //x = x - ((x - (nx - rad)) * (x >= (nx - rad)));
        //y = y - ((y - (ny - rad)) * (y >= (ny - rad)));
        //z = z - ((z - (nz - rad)) * (z >= (nz - rad)));
        	      
         o = orie[z * ny * nx + y * nx + x];
         if (o <= 10)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
					l = k * nx * ny + j * nx + i;
           				//if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde1[l] == 0.0 && fielde2[l] == 0.0 && fielde3[l] == 0.0 && fielde4[l] == 0.0 && fielde5[l] == 0.0 && fielde6[l] == 0.0 && fielde7[l] == 0.0 && fielde8[l] == 0.0 && fielde9[l] == 0.0 && fielde10[l] == 0.0 && fielde11[l] == 0.0 && fielde12[l] == 0.0)
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde1[l] == 0.0)
           				{
           					fielde1[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 10 && o <= 20)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde2[l] == 0.0)
           				{
           					fielde2[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 20 && o <= 30)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde3[l] == 0.0)
           				{
           					fielde3[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 30 && o <= 40)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde4[l] == 0.0)
           				{
           					fielde4[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 40 && o <= 50)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde5[l] == 0.0)
           				{
           					fielde5[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 50 && o <= 60)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde6[l] == 0.0)
           				{
           					fielde6[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 60 && o <= 70)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde7[l] == 0.0)
           				{
           					fielde7[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 70 && o <= 80)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde8[l] == 0.0)
           				{
           					fielde8[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 80 && o <= 90)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde9[l] == 0.0)
           				{
           					fielde9[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 90 && o <= 100)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde10[l] == 0.0)
           				{
           					fielde10[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 100 && o <= 110)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde11[l] == 0.0)
           				{
           					fielde11[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
         else if (o > 110 && o <= 125)
         {      
                 for (i = ix1 + rad; i < ix4 - rad; i++)
                 {
                 	for (j = iy1 + rad; j < iy4 - rad; j++)
           	 	{
           			for (k = iz1 + rad; k < iz4 - rad; k++)
           		  	{
           				if (sqrt((i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z)) < rad && fielde12[l] == 0.0)
           				{
           					fielde12[k * nx * ny + j * nx + i] = 1.0;
           				}				
           		  	}
           	    	}
                 }
         }
  }
}
*/

void initial_ac (double *fieldep, double *fieldeq)
{
  double thrp, thrq;
  double rad = 10.0;
  int i, j, k;
  int x1, y1, z1, l;
  thrp = rand() / (double)(RAND_MAX);
  thrq = rand() / (double)(RAND_MAX);
  if (thrp < 0.01)
  {
    x1 = rad + (rand() / (double)(RAND_MAX)) * (nx - 2 * rad - 1);
    y1 = rad + (rand() / (double)(RAND_MAX)) * (ny - 2 * rad - 1);
    z1 = rad + (rand() / (double)(RAND_MAX)) * (nz - 2 * rad - 1);
    for (i = x1 - rad; i < x1 + rad; i++)
    {
     	for (j = y1 - rad; j < y1 + rad; j++)
        {
    	    for (k = z1 - rad; k < z1 + rad; k++)
      	    {
      		    if (sqrt((i - x1) * (i - x1) + (j - y1) * (j - y1) + (k - z1) * (k - z1)) < rad)
                {
                    double s = 0.0;
                    for (l = 0; l < nac; l++)
                    {
                        s += ac[l].fieldE[k * nx * ny + j * nx + i] * ac[l].fieldE[k * nx * ny + j * nx + i]; 
                    }
                    if (s < 0.1)
                    {
                        fieldep[k * nx * ny + j * nx + i] = 1.0;
                    }
                }
                
   		    }				
     	}
   	 }     
  }
  if (thrq < 0.01)
  {
    x1 = rad + (rand() / (double)(RAND_MAX)) * (nx - 2 * rad - 1);
    y1 = rad + (rand() / (double)(RAND_MAX)) * (ny - 2 * rad - 1);
    z1 = rad + (rand() / (double)(RAND_MAX)) * (nz - 2 * rad - 1);
    for (i = x1 - rad; i < x1 + rad; i++)
    {
     	for (j = y1 - rad; j < y1 + rad; j++)
        {
    	    for (k = z1 - rad; k < z1 + rad; k++)
      	    {
      		    if (sqrt((i - x1) * (i - x1) + (j - y1) * (j - y1) + (k - z1) * (k - z1)) < rad)
                {
                    double s = 0.0;
                    for (l = 0; l < nac; l++)
                    {
                        s += ac[l].fieldE[k * nx * ny + j * nx + i] * ac[l].fieldE[k * nx * ny + j * nx + i]; 
                    }
                    if (s < 0.1)
                    {
                        fieldeq[k * nx * ny + j * nx + i] = 1.0;
                    }
                }
                
   		    }				
     	}
   	 }     
  }

}
#if 0
void initial_ac (double *fieldep, double *fieldeq)
{
  double thrp, thrq;
  double rad = 10.0;
  int i, j, k;
  int x1, y1, z1, x2, y2, z2, l;
  thrp = rand() / (double)(RAND_MAX);
  thrq = rand() / (double)(RAND_MAX);
  int o; 

  if (thrp < 0.8)
  {
      x1 = (int)((rand() / (double)(RAND_MAX)) * nx);
      y1 = (int)((rand() / (double)(RAND_MAX)) * ny);
      z1 = (int)((rand() / (double)(RAND_MAX)) * nz);
      
      for (i = ix1 + rad; i < ix4 - rad; i++)
      {
       	for (j = iy1 + rad; j < iy4 - rad; j++)
         {
      	    for (k = iz1 + rad; k < iz4 - rad; k++)
        	{
	  	        l = k * nx * ny + j * nx + i;
      		    if (sqrt((i - x1) * (i - x1) + (j - y1) * (j - y1) + (k - z1) * (k - z1)) < rad && ac[0].fieldE[l] == 0.0 && ac[1].fieldE[l] == 0.0 && ac[2].fieldE[l] == 0.0 && ac[3].fieldE[l] == 0.0 && ac[4].fieldE[l] == 0.0 && ac[5].fieldE[l] == 0.0)
      		    {
      			    fieldep[k * nx * ny + j * nx + i] = 1.0;
      		    }				
        	}
      	 }   
      }
  }
  if (thrq < 0.8)
  {
      x2 = (int)((rand() / (double)(RAND_MAX)) * nx);
      y2 = (int)((rand() / (double)(RAND_MAX)) * ny);
      z2 = (int)((rand() / (double)(RAND_MAX)) * nz);
      
      for (i = ix1 + rad; i < ix4 - rad; i++)
      {
       	for (j = iy1 + rad; j < iy4 - rad; j++)
         {
      	    for (k = iz1 + rad; k < iz4 - rad; k++)
        	{
	  	        l = k * nx * ny + j * nx + i;
      		    if (sqrt((i - x2) * (i - x2) + (j - y2) * (j - y2) + (k - z2) * (k - z2)) < rad && ac[0].fieldE[l] == 0.0 && ac[1].fieldE[l] == 0.0 && ac[2].fieldE[l] == 0.0 && ac[3].fieldE[l] == 0.0 && ac[4].fieldE[l] == 0.0 && ac[5].fieldE[l] == 0.0)
      		    {
      			    fieldeq[k * nx * ny + j * nx + i] = 1.0;
      		    }				
        	}
      	 }   
      }
  }
}
#endif

void initial_ch (double *fielde1, double *fielde2, double *fielde3, double *fielde4)
{
  int i, j, k;
  double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12;


  for (k = iz1; k < iz4; k++)
  {
     for (j = iy1; j < iy4; j++)
     {
       for (i = ix1; i < ix4; i++)
       {
           e1 = fielde1[k * ny * nx + j * nx + i];
           e2 = fielde2[k * ny * nx + j * nx + i];
           e3 = fielde3[k * ny * nx + j * nx + i];
           e4 = fielde4[k * ny * nx + j * nx + i];
	   if (e1 == 1.0 || e2 == 1.0 || e3 == 1.0 || e4 == 1.0)
	   {
	     ch[0].fieldCI[k * ny * nx + j * nx + i] = 0.1057;
	     ch[1].fieldCI[k * ny * nx + j * nx + i] = 0.0223;
	   }
	   else
	   {
	     ch[0].fieldCI[k * ny * nx + j * nx + i] = 0.102;
	     ch[1].fieldCI[k * ny * nx + j * nx + i] = 0.036;
	   }
       }
     }
  }
}

void
init_field() {
  int n;
  char fname[1024];
  if (restart == 1) {
    for (n = 0; n < nac; n++) {
      sprintf (fname, "%sfield_iter%d_ac%d_%d%d%d.txt", "./restart/", restart_iter, n, cart_id[0], cart_id[1], cart_id[2]);
      read_field (fname, ac[n].fieldE);
      hipMemcpy (fieldE + n * offset, ac[n].fieldE, nx * ny * nz * sizeof (double), hipMemcpyHostToDevice);
    }
  } else {
    //init_field_check ();
    if (nch == 0) {
      for (n = 0; n < nac; n++) {
	printf("nch is zone!\n");
      }
    }
    if (nac == 0) {
      for (n = 0; n < nch; n++) {
	printf("nac is zone!\n");
      }
    }
    if (nac > 0 && nch > 0) {
      time_t t;
      //int p, q;
      srand((unsigned)time(&t) * (myrank + 1));
      initial_zero (ac[0].fieldE, ac[1].fieldE, ac[2].fieldE, ac[3].fieldE, ac[4].fieldE, ac[5].fieldE);
      //init_field_origin();
      P = (int)((rand() / (double)(RAND_MAX)) * 6);
      Q = (int)((rand() / (double)(RAND_MAX)) * 6);
      printf ("myrank = %d p=%d q=%d\n", myrank, P, Q);
      //for (int m = 0; m < 8; m++)
      //{
      //  initial_ac (ac[p].fieldE, ac[q].fieldE);
      //}
//      initial_ch (ac[0].fieldE, ac[1].fieldE, ac[2].fieldE, ac[3].fieldE);
//     for (int m = 0; m < nac; m++)
//     {
//        hipMemcpy (fieldE + m * offset, ac[m].fieldE, nx * ny * nz * sizeof (double), hipMemcpyHostToDevice);
//     }
//     for (int m = 0; m < nch; m++)
//     {
//        hipMemcpy (fieldCI + m * offset, ch[m].fieldCI, nx * ny * nz * sizeof (double), hipMemcpyHostToDevice);
//     }
    }
  }
}

