#include "CNN.h"

CNN::CNN(std::string graph_fn, std::string checkpointPath, int rank)
{
  auto options = tensorflow::SessionOptions();
  detect_result = (float *)malloc(sizeof(float) * 1000 * 17);
  // options.config.mutable_gpu_options()->set_visible_device_list(std::to_string(rank));
  // options.config.mutable_gpu_options()->set_visible_device_list("");
  // options.config.mutable_gpu_options()->clear_visible_device_list();

  auto *device_count = options.config.mutable_device_count();
  device_count->insert({"CPU", 1});
  device_count->insert({"GPU", 0});

  Status status = NewSession(options, &session);
  Status status_load = ReadBinaryProto(Env::Default(), graph_fn, &graphdef); //从meta文件中读取图模型;
  if (!status_load.ok())
  {
    std::cout << "ERROR: Loading model failed..." << graph_fn << std::endl;
    std::cout << status_load.ToString() << "\n";
  }

  Status status_create = session->Create(graphdef.graph_def()); //将模型导入会话Session中;
  if (!status_create.ok())
  {
    std::cout << "ERROR: Creating graph in session failed..." << status_create.ToString() << std::endl;
  }
  // cout << " Session successfully created.Load model successfully!" << endl;

  // 读入预先训练好的模型的权重
  Tensor checkpointPathTensor(DT_STRING, TensorShape());
  checkpointPathTensor.scalar<std::string>()() = checkpointPath;
  status = session->Run(
      {
          {graphdef.saver_def().filename_tensor_name(), checkpointPathTensor},
      },
      {}, {graphdef.saver_def().restore_op_name()}, nullptr);

  if (!status.ok())
  {
    throw runtime_error("Error loading checkpoint from " + checkpointPath + ": " + status.ToString());
  }
  // cout << " Load weights successfully!" << endl;
}

void CNN::eval(std::vector<std::pair<string, Tensor>> inputs, vector<tensorflow::Tensor> &outputs, char *save_path)
{
  Status status_run = session->Run(inputs, {output_node_name}, {}, &outputs);
  if (!status_run.ok())
  {
    std::cout << "ERROR: RUN failed..." << std::endl;
    std::cout << status_run.ToString() << "\n";
  }
  // cout << "---------------------推断完成--------------------" << endl;
  // cout << save_path << endl;
  // printf("保存路径:%s\n", save_path);

  // 获取输出tensor值并保存成文本
  Tensor t = outputs[0];            // Fetch the first tensor
  auto tmap = t.tensor<float, 3>(); // Tensor Shape: [batch_size, target_class_num]

  float *topk_result = (float *)malloc(1000 * 17 * sizeof(float));

  float conf, thresshold = 0.5;
  int count = 0;
  for (int i = 0; i < 1000; i++)
  {
    conf = tmap(0, i, 15);
    if (conf < thresshold)
    {
      continue;
    }
    else
    {
      // for (int j = 0; j < 16; j++)
      // {
      //   printf("%f ", tmap(0, i, j));
      // }

      // printf("\n");

      // 置信度
      topk_result[i * 17 + 0] = tmap(0, i, 15);
      // center
      topk_result[i * 17 + 1] = tmap(0, i, 0)*4;
      topk_result[i * 17 + 2] = tmap(0, i, 1)*4;
      topk_result[i * 17 + 3] = tmap(0, i, 2)*4;
      // whd
      topk_result[i * 17 + 4] = tmap(0, i, 3)*4;
      topk_result[i * 17 + 5] = tmap(0, i, 4)*4;
      topk_result[i * 17 + 6] = tmap(0, i, 5)*4;
      // m
      topk_result[i * 17 + 7] = tmap(0, i, 6);
      topk_result[i * 17 + 8] = tmap(0, i, 7);
      topk_result[i * 17 + 9] = tmap(0, i, 8);
      topk_result[i * 17 + 10] = tmap(0, i, 9);
      topk_result[i * 17 + 11] = tmap(0, i, 10);
      topk_result[i * 17 + 12] = tmap(0, i, 11);
      topk_result[i * 17 + 13] = tmap(0, i, 12);
      topk_result[i * 17 + 14] = tmap(0, i, 13);
      topk_result[i * 17 + 15] = tmap(0, i, 14);
      // class
      topk_result[i * 17 + 16] = 0.0;
    }
    count++;
  }

  // printf("%s\n", full_name);

  FILE *fp_out = fopen(save_path, "w"); //打开输出文件
  if (fp_out == NULL)
  {
    printf("Fail to save\n");
    exit(0);
  }

  for (int i = 0; i < count; i++)
  {
    // printf("%f\n", (output_buffer + ((buffer_i + 1) % 2) * buffer_block_size)[i]);
    // fprintf(fp_out, "%f\n", (output_buffer + ((buffer_i + 1) % 2) * buffer_block_size)[i]);

    // 分别是cls, x, y, z, w, h, d, m*9, conf
    fprintf(fp_out, "[%2.1f, %f, %f, %f, %f, %f, %f, %f, %f ,%f, %f, %f, %f, %f, %f, %f, %f]\n",
            0.0, topk_result[i * 17 + 1], topk_result[i * 17 + 2], topk_result[i * 17 + 3], topk_result[i * 17 + 4], topk_result[i * 17 + 5], topk_result[i * 17 + 6],
            topk_result[i * 17 + 7], topk_result[i * 17 + 8], topk_result[i * 17 + 9], topk_result[i * 17 + 10], topk_result[i * 17 + 11], topk_result[i * 17 + 12], topk_result[i * 17 + 13], topk_result[i * 17 + 14], topk_result[i * 17 + 15],
            topk_result[i * 17 + 0]); // res_4_conv2
  }
  fclose(fp_out);

  // cout << "---------------------后处理完成--------------------" << endl;
}

void CNN::read_dat(std::string dat_path, Tensor *input_tensor)
{
  char str[255];

  float *p = input_tensor->flat<float>().data();

  if ((fp = fopen(dat_path.c_str(), "rt")) == NULL)
  {
    puts("Fail to open file!");
    exit(0);
  }

  for (int i = 0; fgets(str, 255, fp) != NULL; i++)
  {
    p[i] = atof(str);
    if (p[i] > 0.1)
      p[i] = 1;
    else
      p[i] = 0;
  }

  fclose(fp);
}

void CNN::read_dat2(std::string dat_path, float *input_dat)
{
  char str[255];

  if ((fp = fopen(dat_path.c_str(), "rt")) == NULL)
  {
    puts("Fail to open file!");
    exit(0);
  }

  for (int i = 0; fgets(str, 255, fp) != NULL; i++)
  {
    input_dat[i] = atof(str);
    if (input_dat[i] > 0.1)
      input_dat[i] = 1;
    else
      input_dat[i] = 0;
  }

  fclose(fp);
}
