ROCBLAS_DIR = /opt/rocm/rocblas

ROCBLAS_INC = $(ROCBLAS_DIR)/include

ROCBLAS_LIB = $(ROCBLAS_DIR)/lib

HIP_INC     = /opt/rocm/hip/include

HIP_LIB     = /opt/rocm/hip/lib

TENSORFLOW_DIR = /public/home/cnicfengzc/fengzc/tensorflow-rocm/tensorflow-upstream-r1.14-rocm-xuan

TENSORFLOW_INCLUDES = -I $(TENSORFLOW_DIR) -I $(TENSORFLOW_DIR)/bazel-genfiles -I $(TENSORFLOW_DIR)/bazel-bin/tensorflow

OTHER_INCLUDE = -I/public/home/cnicfengzc/fengzc/tensorflow-rocm/tensorflow-upstream-r1.14-rocm-xuan/bazel-genfiles/external/local_config_rocm/rocm/rocm/include/ -I/opt/hpc/software/mpi/hpcx/v2.4.1/gcc-7.3.1/include

PROTOBUF_INCLUDE = -I/public/home/cnicfengzc/fengzc/project_tf_cpp/protobuf/include

EIGEN_INCLUDE = -I/public/home/cnicfengzc/fengzc/usr/local/include/eigen3

ABSAL_INCLUDE = -I/public/home/cnicfengzc/fengzc/usr/local/include/absl

TF_INCLUDES = $(TENSORFLOW_INCLUDES) $(OTHER_INCLUDE) $(EIGEN_INCLUDE) $(ABSAL_INCLUDE) $(PROTOBUF_INCLUDE)

TF_LIBS = -L "./lib" -ltensorflow_cc -ltensorflow_framework

PYTHON_INCLUDE = -I/usr/local/python382/include/python3.8

NUMPY_INCLUDE_DIRS = -I/public/home/cnicfengzc/.local/lib/python3.6/site-packages/numpy/core/include/numpy

PY_INCLUDES = $(NUMPY_INCLUDE_DIRS) $(PYTHON_INCLUDE)

PY_LIBS = -L "/work1/cnicfengzc/project_phase_field/eta12_c2_dgemm_elastic_GX_noelastic_AI/lib" -lpython3.6m

DIR_SRC     = ./src

DIR_OBJ     = ./obj

DIR_BIN     = ./bin

SOURCES_H   = $(wildcard ${DIR_SRC}/*.h)

SOURCES     = $(wildcard ${DIR_SRC}/*.c)

OBJECTS     = $(patsubst %.c,${DIR_OBJ}/%.o,$(notdir ${SOURCES}))

PROGS       = ${DIR_BIN}/ScLETD

#CC_H        = OMPI_CXX=hipcc /opt/hpc/software/mpi/hpcx/v2.4.1/gcc-7.3.1/bin/mpicxx
CC_H        = OMPI_CXX=hipcc mpic++ --amdgpu-target=gfx906

LD          = /opt/hpc/software/mpi/hpcx/v2.7.4/gcc-7.3.1/bin/mpicxx

CFLAGS_H    = -O3 -I$(ROCBLAS_INC) -I$(HIP_INC) 

CFLAGS_D    = -O3 -I$(ROCBLAS_INC) -I$(HIP_INC) 

LDFLAGS     = -lm -L$(ROCBLAS_LIB) -lrocblas  -L$(HIP_LIB) -lhip_hcc 

${PROGS} : ${OBJECTS} ./obj/CNN.o ./obj/ScLETD.o
	$(CC_H) $(OBJECTS) ./obj/CNN.o ./obj/ScLETD.o -o $@ $(LDFLAGS) $(PY_INCLUDES) $(PY_LIBS)  $(TF_INCLUDES) $(TF_LIBS) -w
#$(LD)  $(OBJECTS) ./obj/CNN.o -o $@ $(LDFLAGS) -w $(TF_INCLUDES) $(TF_LIBS)

${DIR_OBJ}/%.o : ${DIR_SRC}/%.c
	$(CC_H)  $(CFLAGS_H) -c  $< -o $@ -w 

./obj/ScLETD.o : ./src/ScLETD.cpp ./src/ScLETD.h
	$(LD) -c ./src/ScLETD.cpp -o ./obj/ScLETD.o $(PY_INCLUDES) $(PY_LIBS) $(TF_INCLUDES) $(TF_LIBS)  

./obj/CNN.o : ./src/CNN.cpp ./src/CNN.h
	$(LD) -c ./src/CNN.cpp -o ./obj/CNN.o  $(PY_INCLUDES) $(PY_LIBS)  $(TF_INCLUDES) $(TF_LIBS)

run:
	mpirun -n 1 ./bin/ScLETD 64

.PHONY:clean
clean:
	rm -f ./obj/*
