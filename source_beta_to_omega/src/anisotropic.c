#include <stdio.h>
#include "mpi.h"
#include "ScLETD.h"
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
// f1
void ac_calc_F1(int n)
{
	int m;
	int i, j, k;
	int n_left, n_right, n_top, n_bottom, n_front, n_back;
	double f_left, f_right, f_top, f_bottom, f_front, f_back;
	double tmp1, tmp2;
  double u0, u1, u2, u3, un;
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				SWTICH_AC_FIELD1(ac[n].f1[k * nx * ny + j * nx + i]);
			}
		}
	}
/*
//grain bundary
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix1 + nghost; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += 10.0 * ac[n].fieldE[k * nx * ny + j * nx + i];
			}
		}
	}
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix4 - 2; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += 10.0 * ac[n].fieldE[k * nx * ny + j * nx + i];
			}
		}
	}
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy1 + nghost; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += 10.0 * ac[n].fieldE[k * nx * ny + j * nx + i];
			}
		}
	}
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy4 - 2; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += 10.0 * ac[n].fieldE[k * nx * ny + j * nx + i];
			}
		}
	}
	for (k = iz1; k < iz1 + nghost; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += 10.0 * ac[n].fieldE[k * nx * ny + j * nx + i];
			}
		}
	}
	for (k = iz4 - 2; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += 10.0 * ac[n].fieldE[k * nx * ny + j * nx + i];
			}
		}
	}*/
//grain bundary
	if (left >= 0)
	{
		i = ix1;
		for (k = iz1; k < iz4; k++)
		{
			for (j = iy1; j < iy4; j++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += ac[n].LE * epn2 * ac[n].fieldEu_left[k * ny + j] / hx / hx;
			}
		}
	}

	if (right >= 0)
	{
		i = ix4 - 1;
		for (k = iz1; k < iz4; k++)
		{
			for (j = iy1; j < iy4; j++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += ac[n].LE * epn2 * ac[n].fieldEu_right[k * ny + j] / hx / hx;
			}
		}
	}

	if (top >= 0)
	{
		j = iy1;
		for (k = iz1; k < iz4; k++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += ac[n].LE * epn2 * ac[n].fieldEu_top[k * nx + i] / hy / hy;
			}
		}
	}

	if (bottom >= 0)
	{
		j = iy4 - 1;
		for (k = iz1; k < iz4; k++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += ac[n].LE * epn2 * ac[n].fieldEu_bottom[k * nx + i] / hy / hy;
			}
		}
	}

	if (front >= 0)
	{
		k = iz1;
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += ac[n].LE * epn2 * ac[n].fieldEu_front[j * nx + i] / hz / hz;
			}
		}
	}

	if (back >= 0)
	{
		k = iz4 - 1;
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ac[n].f1[k * nx * ny + j * nx + i] += ac[n].LE * epn2 * ac[n].fieldEu_back[j * nx + i] / hz / hz;
			}
		}
	}
}

// fieldEt = divergence(D_lambda*gradient(eta))
void ac_calc_F2(n)
{
	int i, j, k;
	int n_left, n_right, n_top, n_bottom, n_front, n_back;
	double a_left, a_right, a_top, a_bottom, a_front, a_back, a_middle;
	double m_left, m_right, m_top, m_bottom, m_front, m_back, m_middle;
	// D_lambda * gradient(eta)
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
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
					a_front = ac[n].fieldE[n_front * nx * ny + j * nx + i];
				}
				else
				{
					a_front = 0;
				}
				if (n_back > -1)
				{
					a_back = ac[n].fieldE[n_back * nx * ny + j * nx + i];
				}
				else
				{
					a_back = 0;
				}
				if (n_top > -1)
				{
					a_top = ac[n].fieldE[k * nx * ny + n_top * nx + i];
				}
				else
				{
					a_top = 0;
				}
				if (n_bottom > -1)
				{
					a_bottom = ac[n].fieldE[k * nx * ny + n_bottom * nx + i];
				}
				else
				{
					a_bottom = 0;
				}
				if (n_right > -1)
				{
					a_right = ac[n].fieldE[k * nx * ny + j * nx + n_right];
				}
				else
				{
					a_right = 0;
				}
				if (n_left > -1)
				{
					a_left = ac[n].fieldE[k * nx * ny + j * nx + n_left];
				}
				else
				{
					a_left = 0;
				}
				// gradient(eta)
				double gradx = (a_right - a_left) / (2.0 * hx);
				double grady = (a_top - a_bottom) / (2.0 * hy);
				double gradz = (a_front - a_back) / (2.0 * hz);
				ac[n].gradx[k * nx * ny + j * nx + i] = ac[n].lambda[0][0] * gradx + ac[n].lambda[0][1] * grady + ac[n].lambda[0][2] * gradz;
				ac[n].grady[k * nx * ny + j * nx + i] = ac[n].lambda[1][0] * gradx + ac[n].lambda[1][1] * grady + ac[n].lambda[1][2] * gradz;
				ac[n].gradz[k * nx * ny + j * nx + i] = ac[n].lambda[2][0] * gradx + ac[n].lambda[2][1] * grady + ac[n].lambda[2][2] * gradz;
			}
		}
	}

	// divergence(D_lambda*gradient_eta)
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
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
					a_front = ac[n].gradz[n_front * nx * ny + j * nx + i];
				}
				else
				{
					a_front = 0;
				}
				if (n_back > -1)
				{
					a_back = ac[n].gradz[n_back * nx * ny + j * nx + i];
				}
				else
				{
					a_back = 0;
				}
				if (n_top > -1)
				{
					a_top = ac[n].grady[k * nx * ny + n_top * nx + i];
				}
				else
				{
					a_top = 0;
				}
				if (n_bottom > -1)
				{
					a_bottom = ac[n].grady[k * nx * ny + n_bottom * nx + i];
				}
				else
				{
					a_bottom = 0;
				}
				if (n_right > -1)
				{
					a_right = ac[n].gradx[k * nx * ny + j * nx + n_right];
				}
				else
				{
					a_right = 0;
				}
				if (n_left > -1)
				{
					a_left = ac[n].gradx[k * nx * ny + j * nx + n_left];
				}
				else
				{
					a_left = 0;
				}
				double tmp = 0.0;
				tmp += (a_right - a_left) / 2.0 / hx;
				tmp += (a_top - a_bottom) / 2.0 / hy;
				tmp += (a_front - a_back) / 2.0 / hz;
				ac[n].f2[k * nx * ny + j * nx + i] = tmp;
			}
		}
	}
}

// divergence((D-lambda)gradient(eta))
void anisotropic_input()
{
	int i, j, n;
	double D[3][3];
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

	fp = fopen("anisotropic_input.txt", "r");
	if (fp == NULL)
	{
		printf("open anisotropic_input.txt Fail!");
		exit(EXIT_FAILURE);
	}

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
