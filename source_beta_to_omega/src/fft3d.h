#ifndef _FFT3D_H_
#define _FFT3D_H_
#include "fft3d_wrap.h"
//#include </public/sourcecode/lammps-30Oct14/src/fft3d_wrap.h>

#define FFT_MEMALIGN 64
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

int xlo, xhi, ylo, yhi, zlo, zhi; // initial partition of grid
int nfft, nfftx, nffty, nfftz;    // size of grid
int fftsize;                      // FFT buffer size returned by FFT setup
void *fft;

FFT_SCALAR *fft_in;
FFT_SCALAR *fft_out;

void fft_setup();
void fft_forward(double *in, double *out_re, double *out_im);
void fft_backward(double *in_re, double *in_im, double *out);
void fft_finish();
void fft_check(double *in, double *out);
void fft_test();

#endif
