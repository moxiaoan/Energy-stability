#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "ScLETD.h"
#include "mkl.h"

/*
la(i,j,k)=(a(i-1,j,k)+a(i+1,j,k)-2*a(i,j,k))/hx^2
        +(a(i,j-1,k)+a(i,j+1,k)-2*a(i,j,k))/hy^2
        +(a(i,j,k-1)+a(i,j,k+1)-2*a(i,j,k))/hz^2;
*/
void laplace(double *la, double *a)
{
	int i, j, k;
	int n_left, n_right, n_top, n_bottom, n_front, n_back;
	double a_left, a_right, a_top, a_bottom, a_front, a_back;
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
					a_front = a[n_front * nx * ny + j * nx + i];
				}
				else
				{
					a_front = 0;
				}
				if (n_back > -1)
				{
					a_back = a[n_back * nx * ny + j * nx + i];
				}
				else
				{
					a_back = 0;
				}
				if (n_top > -1)
				{
					a_top = a[k * nx * ny + n_top * nx + i];
				}
				else
				{
					a_top = 0;
				}
				if (n_bottom > -1)
				{
					a_bottom = a[k * nx * ny + n_bottom * nx + i];
				}
				else
				{
					a_bottom = 0;
				}
				if (n_right > -1)
				{
					a_right = a[k * nx * ny + j * nx + n_right];
				}
				else
				{
					a_right = 0;
				}
				if (n_left > -1)
				{
					a_left = a[k * nx * ny + j * nx + n_left];
				}
				else
				{
					a_left = 0;
				}
				double tmp = 0.0;
				tmp += (a_left + a_right - 2.0 * a[k * nx * ny + j * nx + i]) / hx / hx;
				tmp += (a_top + a_bottom - 2.0 * a[k * nx * ny + j * nx + i]) / hy / hy;
				tmp += (a_front + a_back - 2.0 * a[k * nx * ny + j * nx + i]) / hz / hz;
				la[k * nx * ny + j * nx + i] = tmp;
			}
		}
	}
}

// get f of nth functin
void ch_calc_F0(int n)
{
	int i, j, k, m;
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				SWITCH_CH_FIELD1(ch[n].f0[k * nx * ny + j * nx + i]);
			}
		}
	}
}

//f1 = -C * laplace(-f + KU + Wu) - C * Wmu
void ch_calc_F1(int n)
{
	int i, j, k;
	// ft = -f + KU
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ch[n].ft[k * nx * ny + j * nx + i] = -ch[n].f0[k * nx * ny + j * nx + i] + ch[n].KCI * ch[n].fieldCI[k * nx * ny + j * nx + i];
			}
		}
	}
	if (left >= 0)
	{
		i = ix1;
		for (k = iz1; k < iz4; k++)
		{
			for (j = iy1; j < iy4; j++)
			{
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_left[k * ny + j] / hx / hx;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_right[k * ny + j] / hx / hx;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_top[k * nx + i] / hy / hy;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_bottom[k * nx + i] / hy / hy;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_front[j * nx + i] / hz / hz;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_back[j * nx + i] / hz / hz;
			}
		}
	}
	// f1 = laplace(ft)
	laplace(ch[n].f1, ch[n].ft);
	// f1 = -C * f1
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ch[n].f1[k * nx * ny + j * nx + i] = -ch[n].C * ch[n].f1[k * nx * ny + j * nx + i];
			}
		}
	}
	// f1 -= C * Wmu
	if (left >= 0)
	{
		i = ix1;
		for (k = iz1; k < iz4; k++)
		{
			for (j = iy1; j < iy4; j++)
			{
				ch[n].f1[k * nx * ny + j * nx + i] -= ch[n].C * ch[n].fieldCImu_left[k * ny + j] / hx / hx;
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
				ch[n].f1[k * nx * ny + j * nx + i] -= ch[n].C * ch[n].fieldCImu_right[k * ny + j] / hx / hx;
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
				ch[n].f1[k * nx * ny + j * nx + i] -= ch[n].C * ch[n].fieldCImu_top[k * nx + i] / hy / hy;
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
				ch[n].f1[k * nx * ny + j * nx + i] -= ch[n].C * ch[n].fieldCImu_bottom[k * nx + i] / hy / hy;
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
				ch[n].f1[k * nx * ny + j * nx + i] -= ch[n].C * ch[n].fieldCImu_front[j * nx + i] / hz / hz;
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
				ch[n].f1[k * nx * ny + j * nx + i] -= ch[n].C * ch[n].fieldCImu_back[j * nx + i] / hz / hz;
			}
		}
	}
}

// f2 = divergence(M.*gradient(ft))
void ch_divergence_M_gradient(int n)
{
	int i, j, k;
	int n_left, n_right, n_top, n_bottom, n_front, n_back;
	double a_left, a_right, a_top, a_bottom, a_front, a_back, a_middle;
	double m_left, m_right, m_top, m_bottom, m_front, m_back, m_middle;
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
					a_front = ch[n].ft[n_front * nx * ny + j * nx + i];
				}
				else
				{
					a_front = 0;
				}
				if (n_back > -1)
				{
					a_back = ch[n].ft[n_back * nx * ny + j * nx + i];
				}
				else
				{
					a_back = 0;
				}
				if (n_top > -1)
				{
					a_top = ch[n].ft[k * nx * ny + n_top * nx + i];
				}
				else
				{
					a_top = 0;
				}
				if (n_bottom > -1)
				{
					a_bottom = ch[n].ft[k * nx * ny + n_bottom * nx + i];
				}
				else
				{
					a_bottom = 0;
				}
				if (n_right > -1)
				{
					a_right = ch[n].ft[k * nx * ny + j * nx + n_right];
				}
				else
				{
					a_right = 0;
				}
				if (n_left > -1)
				{
					a_left = ch[n].ft[k * nx * ny + j * nx + n_left];
				}
				else
				{
					a_left = 0;
				}

				if (n_front > -1)
				{
					m_front = ch[n].M[n_front * nx * ny + j * nx + i];
				}
				else
				{
					m_front = 0;
				}
				if (n_back > -1)
				{
					m_back = ch[n].M[n_back * nx * ny + j * nx + i];
				}
				else
				{
					m_back = 0;
				}
				if (n_top > -1)
				{
					m_top = ch[n].M[k * nx * ny + n_top * nx + i];
				}
				else
				{
					m_top = 0;
				}
				if (n_bottom > -1)
				{
					m_bottom = ch[n].M[k * nx * ny + n_bottom * nx + i];
				}
				else
				{
					m_bottom = 0;
				}
				if (n_right > -1)
				{
					m_right = ch[n].M[k * nx * ny + j * nx + n_right];
				}
				else
				{
					m_right = 0;
				}
				if (n_left > -1)
				{
					m_left = ch[n].M[k * nx * ny + j * nx + n_left];
				}
				else
				{
					m_left = 0;
				}
				a_middle = ch[n].ft[k * nx * ny + j * nx + i];
				m_middle = ch[n].M[k * nx * ny + j * nx + i];
				double tmp = 0.0;
				tmp += (m_left + m_middle) * (a_left - a_middle) / 2.0 / hx / hx;
				tmp += (m_right + m_middle) * (a_right - a_middle) / 2.0 / hx / hx;
				tmp += (m_top + m_middle) * (a_top - a_middle) / 2.0 / hy / hy;
				tmp += (m_bottom + m_middle) * (a_bottom - a_middle) / 2.0 / hy / hy;
				tmp += (m_front + m_middle) * (a_front - a_middle) / 2.0 / hz / hz;
				tmp += (m_back + m_middle) * (a_back - a_middle) / 2.0 / hz / hz;
				ch[n].f2[k * nx * ny + j * nx + i] = tmp;
			}
		}
	}
}

// f2 = -divergence((M-C)gradient(epn2 * laplace(U) - f + epn2 * Wu))
void ch_calc_F2(int n)
{
	int i, j, k;
	// ft = laplace(U)
	laplace(ch[n].ft, ch[n].fieldCI);
	// ft = epn2 * ft - f
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ch[n].ft[k * nx * ny + j * nx + i] = epn2 * ch[n].ft[k * nx * ny + j * nx + i] - ch[n].f0[k * nx * ny + j * nx + i];
			}
		}
	}
	// ft = ft + epn2 * Wu
	if (left >= 0)
	{
		i = ix1;
		for (k = iz1; k < iz4; k++)
		{
			for (j = iy1; j < iy4; j++)
			{
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_left[k * ny + j] / hx / hx;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_right[k * ny + j] / hx / hx;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_top[k * nx + i] / hy / hy;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_bottom[k * nx + i] / hy / hy;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_front[j * nx + i] / hz / hz;
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
				ch[n].ft[k * nx * ny + j * nx + i] += epn2 * ch[n].fieldCIu_back[j * nx + i] / hz / hz;
			}
		}
	}
	// M = M - C
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ch[n].M[k * nx * ny + j * nx + i] = ch[n].M[k * nx * ny + j * nx + i] - ch[n].C;
			}
		}
	}
	// f2 = divergence(M(gradient(ft)))
	ch_divergence_M_gradient(n);
	// f2 = -f2
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				ch[n].f2[k * nx * ny + j * nx + i] = -ch[n].f2[k * nx * ny + j * nx + i];
			}
		}
	}
}

void ch_calc_FU_VariableMobility(int n, double *f)
{
	int i, j, k;
	// f = ...
	ch_calc_F0(n);
	// f1 = ...
	ch_calc_F1(n);
	// f2 = ...
	ch_calc_F2(n);
	// U1 = f1 + f2
	for (k = iz1; k < iz4; k++)
	{
		for (j = iy1; j < iy4; j++)
		{
			for (i = ix1; i < ix4; i++)
			{
				f[k * nx * ny + j * nx + i] = ch[n].LCI * (ch[n].f1[k * nx * ny + j * nx + i] + ch[n].f2[k * nx * ny + j * nx + i]);
			}
		}
	}
}

void update_M_C()
{
	int i, j, k;
	double cu0, mn0, ni0;

	double QACu, QGCu, D0ACu, D0GCu, QAMn, QGMn, D0AMn, D0GMn, QANi, QGNi, D0ANi, D0GNi;
	double DCuA, DCuG, DMnA, DMnG, DNiA, DNiG;
	double gconst, tempr, RT;
	gconst = 8.314472;
	tempr = 823.0;
	RT = gconst * tempr;
	cu0 = 0.15;
	mn0 = 0.01;
	ni0 = 0.01;
	ch[0].C = -1.0e30;
	ch[1].C = -1.0e30;
	ch[2].C = -1.0e30;

	QACu = 2.44e5;
	QGCu = 2.80e5;

	D0ACu = 4.7e-5;
	D0GCu = 4.3e-5;

	QAMn = 2.63e5;
	QGMn = 2.64e5;

	D0AMn = 1.49e-4;
	D0GMn = 2.78e-5;

	QANi = 2.56e5;
	QGNi = 2.73e5;

	D0ANi = 1.4e-4;
	D0GNi = 1.08e-5;

	DCuA = (D0ACu * exp(-QACu / RT));
	DCuG = (D0GCu * exp(-QGCu / RT)) / DCuA;

	DMnA = (D0AMn * exp(-QAMn / RT)) / DCuA;
	DMnG = (D0GMn * exp(-QGMn / RT)) / DCuA;

	DNiA = (D0ANi * exp(-QANi / RT)) / DCuA;
	DNiG = (D0GNi * exp(-QGNi / RT)) / DCuA;

	DCuA = 1.0;

	for (k = 0; k < nz; k++)
	{
		for (j = 0; j < ny; j++)
		{
			for (i = 0; i < nx; i++)
			{
				ch[0].M[k * ny * nx + j * nx + i] = cu0 * (1.0 - cu0) * (DCuA * (1.0 - ac[0].fieldE[k * ny * nx + j * nx + i]) + DCuG * ac[0].fieldE[k * ny * nx + j * nx + i]);
				ch[1].M[k * ny * nx + j * nx + i] = mn0 * (1.0 - mn0) * (DMnA * (1.0 - ac[0].fieldE[k * ny * nx + j * nx + i]) + DMnG * ac[0].fieldE[k * ny * nx + j * nx + i]);
				ch[2].M[k * ny * nx + j * nx + i] = ni0 * (1.0 - ni0) * (DNiA * (1.0 - ac[0].fieldE[k * ny * nx + j * nx + i]) + DNiG * ac[0].fieldE[k * ny * nx + j * nx + i]);
				ch[0].C = (ch[0].C > ch[0].M[k * ny * nx + j * nx + i]) ? ch[0].C : ch[0].M[k * ny * nx + j * nx + i];
				ch[1].C = (ch[1].C > ch[1].M[k * ny * nx + j * nx + i]) ? ch[1].C : ch[1].M[k * ny * nx + j * nx + i];
				ch[2].C = (ch[2].C > ch[2].M[k * ny * nx + j * nx + i]) ? ch[2].C : ch[2].M[k * ny * nx + j * nx + i];
			}
		}
	}
}

void limit_U()
{
	int i, j, k;
	int n;
	for (n = 0; n < nch; n++)
	{
		for (k = 0; k < nz; k++)
		{
			for (j = 0; j < ny; j++)
			{
				for (i = 0; i < nx; i++)
				{
					if (ch[n].fieldCI[k * ny * nx + j * nx + i] >= 0.9999)
						ch[n].fieldCI[k * ny * nx + j * nx + i] = 0.9999;
					if (ch[n].fieldCI[k * ny * nx + j * nx + i] < 0.00001)
						ch[n].fieldCI[k * ny * nx + j * nx + i] = 0.00001;
				}
			}
		}
	}
}
