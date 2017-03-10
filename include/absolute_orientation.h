#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
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
 void time_synchronization(const vector<PoseDatas> &key_dataset, const vector<PoseDatas> &ground_dataset);
 void compute_q(const vector<PoseDatas> &key_dataset, const vector<PoseDatas> &ground_dataset);
 void compute_s(const vector<Eigen::Matrix<double,3,1>> &rl,const vector<Eigen::Matrix<double,3,1>> &rr);
 void compute_p(const double &s,const Eigen::Quaterniond &q,const Eigen::Matrix<double,3,1> &rl_ave,const Eigen::Matrix<double,3,1> &rr_ave);
 vector<PoseDatas> key_dataset;
 vector<PoseDatas> ground_dataset;
 vector<PoseDatas> some_ground_dataset;
 vector<PoseDatas> temporary_dataset;
 private:
// PoseDatas pose;

 int n; //sum of dataset
 Eigen::Matrix<double,3,1>  rl_ave,rr_ave;
 Eigen::Matrix<double,3,3> M;
 Eigen::Matrix<double,4,4> N;
 Eigen::Quaterniond q;
};
