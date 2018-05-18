#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
#include "matrix.h"
using namespace std;
struct PoseDatas
{
 double time; 
 Eigen::Matrix<double,3,1> position;
};
class Coordinate_Align
{
 public:
 
 void input_datas(const string &keyframe_pose,const string &groundtruth_pose);
 void time_synchronization( vector<PoseDatas> &key_dataset,  vector<PoseDatas> &ground_dataset);
 void compute_q(const vector<PoseDatas> &key_dataset, const vector<PoseDatas> &ground_dataset);
 void compute_s(const vector<Eigen::Matrix<double,3,1>> &rl,const vector<Eigen::Matrix<double,3,1>> &rr);
 void compute_p(const double &s,const Eigen::Quaterniond &q,const Eigen::Matrix<double,3,1> &rl_ave,const Eigen::Matrix<double,3,1> &rr_ave);

 vector<PoseDatas> key_dataset;
 vector<PoseDatas> ground_dataset;
 vector<PoseDatas> some_ground_dataset;
 vector<PoseDatas> temporary_dataset;
 private:
// PoseDatas pose;
 void savePathPlot (vector<Matrix> &poses_gt,vector<Matrix> &poses_result,string file_name) ;
 vector<int32_t> computeRoi (vector<Matrix> &poses_gt,vector<Matrix> &poses_result);
 void plotPathPlot (string dir,vector<int32_t> &roi,int32_t idx);
 int n; //sum of dataset
 Eigen::Matrix<double,3,1>  rl_ave,rr_ave;
 Eigen::Matrix<double,3,3> M;
 Eigen::Matrix<double,4,4> N;
 Eigen::Quaterniond q;

 string file_name;
};
