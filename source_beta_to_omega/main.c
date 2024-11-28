#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif

int main()
{
	int nx = 68;
	int ny = 68;
	int nz = 68;
	double *fieldE;
	double *fieldE1;
	int i, j, k;
	fieldE = (double *)_mm_malloc(sizeof(double) * nx * ny * nz * 8, 256);
	fieldE1 = (double *)_mm_malloc(sizeof(double) * nx * ny * nz, 256);
	char filename[1024];
	FILE *file;
/*	sprintf (filename, "./data/eta3_000007_000000.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + j * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000007_000001.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + j * nx + i - nx];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000007_010000.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + (j-ny) * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000007_010001.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + (j-ny) * nx + i - nx];
        		}
        	}
  	}

	sprintf (filename, "./data/eta3_000007_000100.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + j * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000007_000101.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = 0; j < ny; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + j * nx + i - nx];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000007_010100.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + (j-ny) * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000007_010101.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + (j-ny) * nx + i - nx];
        		}
        	}
  	}*/
	sprintf (filename, "./data/eta3_000002_000000.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + j * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000002_010000.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + j * nx + i - nx];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000002_000100.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + (j-ny) * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000002_010100.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = 0; k < nz; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[k * nx * ny + (j-ny) * nx + i - nx];
        		}
        	}
  	}

	sprintf (filename, "./data/eta3_000002_000001.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + j * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000002_010001.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = 0; j < ny; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + j * nx + i - nx];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000002_000101.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = 0; i < nx; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + (j-ny) * nx + i];
        		}
        	}
  	}
	sprintf (filename, "./data/eta3_000002_010101.dat");
	file = fopen (filename, "r");
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				fscanf (file, "%lf", &fieldE1[k * nx * ny + j * nx + i]);
        		}
        	}
  	}
  	fclose (file);
	for (k = nz; k < nz*2; k++) {
		for (j = ny; j < ny*2; j++) {
			for (i = nx; i < nx*2; i++) {
				fieldE[k * nx*2 * ny*2 + j * nx*2 + i] = fieldE1[(k-nz) * nx * ny + (j-ny) * nx + i - nx];
        		}
        	}
  	}

	sprintf (filename, "./data/change/eta3_000002_010101.dat");
	file = fopen (filename, "w");
	for (k = 0; k < nz*2; k++) {
		for (j = 0; j < ny*2; j++) {
			for (i = 0; i < nx*2; i++) {
				fprintf (file, "%lf\n", fieldE[k * nx*2 * ny*2 + j * nx*2 + i]);
        		}
        	}
  	}
  	fclose (file);
	_mm_free(fieldE);
	_mm_free(fieldE1);
	return 0;
}
