#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"
#include "rocblas.h"
#include <hipfft.h>

hipEvent_t st, ed;
hipEvent_t st2, ed2;
hipDeviceProp_t props;

rocblas_handle handle;
rocblas_int x_m, x_n, x_k;
rocblas_int y_m, y_n, y_k;
rocblas_int z_m, z_n, z_k;
rocblas_int Gx_m, Gx_n, Gx_k;
rocblas_int Gy_m, Gy_n, Gy_k;
rocblas_int Gz_m, Gz_n, Gz_k;