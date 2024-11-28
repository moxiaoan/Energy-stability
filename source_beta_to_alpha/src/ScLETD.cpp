#include "CNN.h" // CNN.h必须放在ScLETD.h的前面否则会出现错误
#include "ScLETD.h"
// int checkpoint_; // 这个必须在main里声明, 不能再ScLETD.h里

int main(int argc, char **argv)
{
  int n;
  const int detect_iter = 10000;
  double tt;
  char filename[1024];
  FILE *file;
  start_mpi(argc, argv);
  sprintf(filename, "input%s.ini", argv[1]);
  file = fopen(filename, "r");
  fscanf(file, "%d,%d,%d", &nx, &ny, &nz);
  fscanf(file, "%d,%d,%d", &procs[0], &procs[1], &procs[2]);
  fscanf(file, "%d", &nghost);
  fscanf(file, "%s", work_dir);
  fscanf(file, "%s", data_dir);
  fscanf(file, "%lf,%lf", &dt, &t_total);
  fscanf(file, "%d", &restart);
  fscanf(file, "%d", &periodic);
  fscanf(file, "%d,   %d", &nac, &nch);
  fscanf(file, "%d", &ELASTIC);
  fscanf(file, "%d", &ANISOTROPIC);
  fscanf(file, "%lf,%lf,%lf", &xmin, &ymin, &zmin);
  fscanf(file, "%lf,%lf,%lf", &xmax, &ymax, &zmax);
  fscanf(file, "%lf,%lf,%lf", &kkx, &kky, &kkz);
  fscanf(file, "%d,%d,%d", &checkpoint_, &nchk, &chk);
  fscanf(file, "%d", &nout);
  fscanf(file, "%d", &Approx);
  fclose(file);

  if (myrank == prank)
  {
    printf("running   on\t%d   processors\n", nprocs);
    printf("restart\t\t%d\n", restart);
    printf("nx,ny,nz\t%d\t%d\t%d\n", nx, ny, nz);
    printf("px,py,pz\t%d\t%d\t%d\n", procs[0], procs[1], procs[2]);
    printf("nghost\t\t%d\n", nghost);
    printf("output   dir\t%s\n", work_dir);
    printf("data   dir\t%s\n", data_dir);
    printf("dt,t_total\t%lf,%lf\n", dt, t_total);
    printf("ELASTIC\t\t%d\n", ELASTIC);
    printf("ANISOTROPIC\t\t%d\n", ANISOTROPIC);
    printf("nac,nch\t\t%d,%d\n", nac, nch);
    printf("xmin,ymin,zmin\t%lf\t%lf\t%lf\n", xmin, ymin, zmin);
    printf("xmax,ymax,zmax\t%lf\t%lf\t%lf\n", xmax, ymax, zmax);
    printf("checkpoint_,nchk,chk\t%d,%d,%d\n", checkpoint_, nchk, chk);
    printf("nout\t\t%d\n", nout);
    printf("Approx\t\t%d\n", Approx);
    printf("detect_iter\t\t%d\n", detect_iter);
  }
  if (nac == 0 && nch == 0)
  {
    printf("nac   =   0   and   nch   =   0\n");
    exit(1);
  }
  if (ELASTIC == 1)
  {
    elastic_input();
  }

  alloc_vars();
  init_vars();

  // MPI_Init(&argv, &argc);
  // int world_size;
  // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string graph_fn = "./CenterNet3D_weight/CenterNet3D_0.0.5.16_v31_7_h2.25to25_20220506-54.16448974609375eval.ckpt.meta";
  const std::string checkpointPath = "./CenterNet3D_weight/CenterNet3D_0.0.5.16_v31_7_h2.25to25_20220506-54.16448974609375eval.ckpt";

  CNN cnn = CNN(graph_fn, checkpointPath, 0);

  tensorflow::Tensor input_tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 128, 128, 128, 1}));
  std::vector<std::pair<string, Tensor>> inputs;
  inputs.emplace_back(std::string("images"), input_tensor);

  float *input_tensor_data = input_tensor.flat<float>().data();

  vector<tensorflow::Tensor> outputs;

  define_mpi_type();
  sw_cart_creat();
  init_para();
    printf("##step%dfieldE, myrank = %d, (%d, %d, %d), on %s\n",iter,myrank,cart_id[0],cart_id[1],cart_id[2],processor_name);
  read_matrices();
  //if (myrank == 0)
  //{
    printf("before read_ori\n");
  //}
  read_ori();
  //if (myrank == 0)
  //{
    printf("after read_ori\n");
  //}
  if (ELASTIC == 1)
  {
    elastic_init();
  }
    printf("before read_chk\n");
  if (checkpoint_ == 1)
  {
    irun = 0;
    chk = 0;
    init_field();
  }
  else if (checkpoint_ == 2)
  {
    read_chk();
    ioutput = 7;
    chk = (chk + 1) % 2;
  }
    printf("after read_chk\n");
    printf("step%dfieldE, myrank = %d, (%d, %d, %d), on %s\n",iter,myrank,cart_id[0],cart_id[1],cart_id[2],processor_name);

  if (ANISOTROPIC == 1)
  {
    anisotropic_input();
  }
  iter = 0;
  tt = 0.0;
  ot1 = 0;
  ot2 = 0;
  out_iter = 0;
  out_tt = 0;
  count_gyq = 0;
  double dt_const = dt;
  while (out_tt < (floor(t_total / dt) / 100))
  {
    out_tt += 1;
    ac_calc_F1(field2_all, iter, detect_iter);
    if (myrank == 0)
    {
        printf("iter = %d, detect_iter = %d\n", iter, detect_iter);
    }
    if (iter % detect_iter == 0)
    {
      for (int m = 0; m < nac; m++)
      {
        char dat_temp_string[128] = "";
        char dat_fullname[256] = "";
        sprintf(dat_temp_string, "%d", iter);
        strcat(dat_fullname, dat_temp_string);
        strcat(dat_fullname, "_");
        sprintf(dat_temp_string, "%d_%d_%d_ds", cart_id[0], cart_id[1], cart_id[2]);
        strcat(dat_fullname, dat_temp_string);
        strcat(dat_fullname, ".txt");

        // 除了保存下采样模拟数据之外还保存对应该步的检测结果
        char filed1_name[1280] = "";
        sprintf(filed1_name, "./downsample_output_%d_%d%d%d/eta%d_", nx, procs[0], procs[1], procs[2], m);
        // char filed1_name[1280] = "./downsample_output_0501_256_222/eta0_";
        strcat(filed1_name, dat_fullname);
        compress_save_output(filed1_name, field2_all + m * 128 * 128 * 128, field2_cp, 128 * 128 * 16);

        char filed1_name_detect[1280] = "";
        sprintf(filed1_name_detect, "./downsample_detect_result_%d_%d%d%d/eta%d_", nx, procs[0], procs[1], procs[2], m);
        strcat(filed1_name_detect, dat_fullname);

        memcpy(input_tensor_data, field2_all + m * 128 * 128 * 128, sizeof(float) * 128 * 128 * 128);

        cnn.eval(inputs, outputs, filed1_name_detect);
      }

      if (myrank == 0)
        printf("eval complete\n");
      MPI_Barrier(MPI_COMM_WORLD);
    }

    out_iter++;
  }
  dealloc_vars();
  close_mpi();

  return 0;
}
