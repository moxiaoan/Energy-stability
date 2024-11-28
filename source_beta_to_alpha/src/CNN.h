#include <iostream>
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/protobuf/meta_graph.pb.h"
#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/core/framework/tensor.h"
#include <fstream>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <mpi.h>
#include <sys/io.h>
#include <time.h>


using namespace std;
using namespace tensorflow;
using namespace tensorflow::ops;
using tensorflow::Status;
using tensorflow::string;
using tensorflow::Tensor;

class CNN
{
public:
  CNN(std::string graph_fn, std::string checkpointPath, int rank);
  void eval(std::vector<std::pair<string, Tensor>> inputs, vector<tensorflow::Tensor> &outputs, char* save_path);
  void read_dat(std::string dat_path, Tensor *input_tensor);
  void read_dat2(std::string dat_path, float *input_dat);
  void getFiles(string path, vector<string> &files);

private:
  Session *session;
  MetaGraphDef graphdef;
  const string input_node_name = "images";
  const string output_node_name = "detections";
  int rank;
  FILE *fp;
  float *detect_result;
  const int input_weight = 128, input_height = 128, input_depth = 128;
};
