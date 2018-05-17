#include "absolute_orientation.h"
#include<iostream>
#include<algorithm>
#include<fstream>
#include <iomanip>
#include <Eigen/Eigenvalues>

#include <stdio.h>
#include <math.h>
#include <limits>
#include <sys/stat.h>
using namespace std;
 void Coordinate_Align::input_datas(const string &keyframe_pose,const string &groundtruth_pose)
 {
     ifstream key_stream,ground_stream;
     string tempStr;
     PoseDatas key_pose,ground_pose;

     key_stream.open(keyframe_pose);
     while(!key_stream.eof())     //acquare keyframe dataset
   {
     getline(key_stream, tempStr);
     if(key_stream.fail())
         break;
     key_stream>>key_pose.time>>key_pose.position(0)>>key_pose.position(1)>>key_pose.position(2);
     key_dataset.push_back(key_pose);
   }
    ground_stream.open(groundtruth_pose);
    while(!ground_stream.eof())   //acquare groundtruth dataset
   {
       getline(ground_stream, tempStr);
        if(ground_stream.fail())
            break;
        ground_stream>>ground_pose.time>>ground_pose.position(0)>>ground_pose.position(1)>>ground_pose.position(2);
        ground_dataset.push_back(ground_pose);
   }
    assert(key_dataset.empty()) ;
    assert(ground_dataset.empty()) ;
  time_synchronization(key_dataset,ground_dataset);

 }
 ///
 /// \brief Coordinate_Align::time_synchronization
 /// \param key_dataset
 /// \param ground_dataset
 ///   acquare groundtruth dataset time-synchronization with keyframe dataset
void  Coordinate_Align::time_synchronization( vector<PoseDatas> &key_dataset,  vector<PoseDatas> &ground_dataset)
{
    //the time of the first keyframe must more than the time of the first ground truth
   while(key_dataset.front().time<=ground_dataset.front().time)
        key_dataset.erase(key_dataset.begin());
 //the time of the last keyframe must less than the time of the last ground truth
   while(key_dataset.back().time>=ground_dataset.back().time)
        key_dataset.erase(key_dataset.end());
  assert(key_dataset.empty());
  int j=0;
  n=key_dataset.size();// duo to my the most time of keyframe dataset more than the most time of groundtruth dataset
  ofstream f;
  f.open("new_groundtruth.txt");
  f << fixed;
  for(int i=0;i<n;i++)
      {
        while(ground_dataset[j].time-key_dataset[i].time<=0)
            {
             temporary_dataset.push_back(ground_dataset[j]);
             j++;

            }
        some_ground_dataset.push_back(temporary_dataset.back());
        f << setprecision(6) <<temporary_dataset.back().time<< " " << setprecision(9) <<temporary_dataset.back().position(0)<< " "
          << temporary_dataset.back().position(1)<< " " << temporary_dataset.back().position(2)<< endl;
      }
  compute_q(key_dataset,some_ground_dataset);
}
void Coordinate_Align::compute_q(const vector<PoseDatas> &key_dataset, const vector<PoseDatas> &ground_dataset)
{
    Eigen::Matrix<double,3,1> rl_err,rr_err;
    vector<Eigen::Matrix<double,3,1>>rl,rr;
    Eigen::Matrix<double,1,3> rr_transpose;
    double Sxx,Sxy,Sxz;
    double Syx,Syy,Syz;
    double Szx,Szy,Szz;
    double max_eigenvalue=0;
    for(int i=0;i<n;i++)
       {
        rl_err=rl_err+key_dataset[i].position;
        rr_err=rr_err+ground_dataset[i].position;
       }
       rl_ave=rl_err/n;
       rr_ave=rr_err/n;
    for(int j=0;j<n;j++)
        {
         rl_err=key_dataset[j].position-rl_ave;
         rr_err=ground_dataset[j].position-rr_ave;
         rl.push_back(rl_err);
         rr.push_back(rr_err);
         }
     for(int i=0; i<n;i++)   //acquare M
        {
         rr_transpose(0)=rr[i](0);
         rr_transpose(1)=rr[i](1);
         rr_transpose(2)=rr[i](2);
        M=M+rl[i]*rr_transpose;
       }
     Sxx=M(0,0);
     Sxy=M(0,1);
     Sxz=M(0,2);
     Syx=M(1,0);
     Syy=M(1,1);
     Syz=M(1,2);
     Szx=M(2,0);
     Szy=M(2,1);
     Szz=M(2,2);

     N(0,0)=Sxx+Syy+Szz;
     N(0,1)=Syz-Szy;
     N(0,2)=Szx-Sxz;
     N(0,3)=Sxy-Syx;
     N(1,0)=Syz-Szy;
     N(1,1)=Sxx-Syy-Szz;
     N(1,2)=Sxy+Syx;
     N(1,3)=Szx+Sxz;
     N(2,0)=Szx-Sxz;
     N(2,1)=Sxy+Syx;
     N(2,2)=-Sxx+Syy-Szz;
     N(2,3)=Syz+Szy;
     N(3,0)=Sxy-Syx;
     N(3,1)=Szx+Sxz;
     N(3,2)=Syz+Szy;
     N(3,3)=-Sxx-Syy+Szz;
     Eigen::EigenSolver<Eigen::Matrix4d> es(N);
     Eigen::Matrix4d D = es.pseudoEigenvalueMatrix();
     Eigen::Matrix4d V = es.pseudoEigenvectors();
     for(int i=0;i<4;i++)
         {
         if(D(i,i)>max_eigenvalue)
         {
         max_eigenvalue= D(i,i);
         q.w()= V(0,i);
         q.vec()(0)=V(1,i);
         q.vec()(1)=V(2,i);
         q.vec()(2)=V(3,i);
         q.normalized();
         }
         }
     cout<<"q= "<<q.w()<<" "<<q.vec()(0)<<" "<<q.vec()(1)<<" "<<q.vec()(2)<<endl;
     compute_s(rl,rr);
     }
void Coordinate_Align::compute_s(const vector<Eigen::Matrix<double,3,1>> &rl,const vector<Eigen::Matrix<double,3,1>> &rr)
{
    double Sr=0,Sl=0,s;
    for(int i=0;i<n;i++)
        {
         Sr=Sr+sqrt(rr[i](0)*rr[i](0)+rr[i](1)*rr[i](1)+rr[i](2)*rr[i](2));
         Sl=Sl+sqrt(rl[i](0)*rl[i](0)+rl[i](1)*rl[i](1)+rl[i](2)*rl[i](2));
        }
    s=sqrt(Sr/Sl);
    cout<<"s= "<<s<<endl;
   compute_p(s,q,rl_ave,rr_ave);

}
void Coordinate_Align::compute_p(const double &s,const Eigen::Quaterniond &q,const Eigen::Matrix<double,3,1> &rl_ave,const Eigen::Matrix<double,3,1> &rr_ave)
{
    Eigen::Matrix<double,3,1> p;
    double error_x,error_y,error_z,error_ave; //estimate value-true value
    //Eigen::Matrix<double,3,1> new_keyframe_pose;
    ofstream f,ff;
    f.open("new_keyframe.txt");
    ff.open("error.txt");
    f << fixed;
    ff << fixed;
    p=rr_ave-s*q.toRotationMatrix()*rl_ave;
    cout<<"p= "<<p<<endl;
    double error_sum =0,total_error_ave;
    for(int i=0;i<n;i++)
      {
       key_dataset[i].position=s*q.toRotationMatrix()*key_dataset[i].position+p;
       error_x=key_dataset[i].position(0)-some_ground_dataset[i].position(0);
       error_y=key_dataset[i].position(1)-some_ground_dataset[i].position(1);
       error_z=key_dataset[i].position(2)-some_ground_dataset[i].position(2);
       error_ave=sqrt(error_x*error_x+error_y*error_y+error_z*error_z);

       f << setprecision(6) <<key_dataset[i].time<< " " << setprecision(9) <<key_dataset[i].position(0)<< " "
         <<key_dataset[i].position(1)<< " " << key_dataset[i].position(2)<< endl;

       ff<< setprecision(6) <<key_dataset[i].time<< " " << setprecision(9)<<error_x<<" "
          <<error_y<<" "<<error_z<<" "<<error_ave<<endl;
       error_sum += error_ave;
      }
     total_error_ave = error_sum/n;
     ff<<"total average error : "<<total_error_ave<<endl;
     cout<<"total_error_ave = "<<total_error_ave<<endl;
  plot_trajectory(key_dataset,some_ground_dataset);
}
void Coordinate_Align::plot_trajectory(const vector<PoseDatas> &key_dataset, const vector<PoseDatas> &ground_dataset)
{
    Matrix tmp_pose (3,1);
    vector<Matrix> poses_result, poses_gt;
   for(auto &keyframe_pose: key_dataset)
   {
       tmp_pose.val[0][0] = keyframe_pose.position(0);
       tmp_pose.val[1][0] = keyframe_pose.position(1);
       tmp_pose.val[2][0] = keyframe_pose.position(2);
       poses_result.push_back(tmp_pose);
   }
   for(auto &ground_pose: ground_dataset)
   {
       tmp_pose.val[0][0] = ground_pose.position(0);
       tmp_pose.val[1][0] = ground_pose.position(1);
       tmp_pose.val[2][0] = ground_pose.position(2);
       poses_gt.push_back(tmp_pose);
   }
  //create file
   char * error_dir  = "errors";
   char * plot_path_dir  =  "plot_path";
   char * plot_error_dir =  "plot_error";
   mkdir(error_dir, 0777);
   mkdir(plot_path_dir, 0777);
   mkdir (plot_error_dir, 0777);
   file_name = "00.txt";

   // check for errors
   if (poses_gt.size()==0 || poses_result.size()!=poses_gt.size())
   {
     cout<<"poses_gt.size()= "<<poses_gt.size()<<endl;
     cout<<"poses_result.size()= "<<poses_result.size()<<endl;
     cout<< "ERROR: Could not read (all) poses of:" << file_name << endl;
     return;
   }
   //save trajectory
   savePathPlot(poses_gt,poses_result, string(plot_path_dir) + "/" + file_name);
   //compute roi
   vector<int32_t> roi = computeRoi(poses_gt,poses_result);
   //plot trajectory
   plotPathPlot(string(plot_path_dir),roi,0);

}
void Coordinate_Align::savePathPlot (vector<Matrix> &poses_gt,vector<Matrix> &poses_result,string file_name)
{
  // open file
  FILE *fp = fopen(file_name.c_str(),"w");
  // save x/z coordinates of all frames to file
  for (int i=0; i<poses_gt.size(); i++)
    fprintf(fp,"%f %f %f %f\n",poses_gt[i].val[0][0],poses_gt[i].val[1][0], poses_result[i].val[0][0],poses_result[i].val[1][0]);
  // close file
  fclose(fp);
}
vector<int32_t> Coordinate_Align::computeRoi (vector<Matrix> &poses_gt,vector<Matrix> &poses_result) {

  float x_min = numeric_limits<int32_t>::max();
  float x_max = numeric_limits<int32_t>::min();
  float z_min = numeric_limits<int32_t>::max();
  float z_max = numeric_limits<int32_t>::min();

  for (vector<Matrix>::iterator it=poses_gt.begin(); it!=poses_gt.end(); it++) {
    float x = it->val[0][1];
    float z = it->val[1][0];
    if (x<x_min) x_min = x; if (x>x_max) x_max = x;
    if (z<z_min) z_min = z; if (z>z_max) z_max = z;
  }

  for (vector<Matrix>::iterator it=poses_result.begin(); it!=poses_result.end(); it++) {
    float x = it->val[0][0];
    float z = it->val[1][0];
    if (x<x_min) x_min = x; if (x>x_max) x_max = x;
    if (z<z_min) z_min = z; if (z>z_max) z_max = z;
  }

  float dx = 1.2*(x_max-x_min);
  float dz = 1.2*(z_max-z_min);
  float mx = 0.6*(x_max+x_min);
  float mz = 0.6*(z_max+z_min);
  float r  = 0.6*max(dx,dz);

  vector<int32_t> roi;
  roi.push_back((int32_t)(mx-r));
  roi.push_back((int32_t)(mx+r));
  roi.push_back((int32_t)(mz-r));
  roi.push_back((int32_t)(mz+r));
  return roi;
}
void Coordinate_Align::plotPathPlot (string dir,vector<int32_t> &roi,int32_t idx) {

  // gnuplot file name
  char command[1024];
  char file_name[256];
  sprintf(file_name,"%02d.gp",idx);
  string full_name = dir + "/" + file_name;
  // create png + eps
  for (int32_t i=0; i<2; i++) {

    // open file
    FILE *fp = fopen(full_name.c_str(),"w");
    // save gnuplot instructions
    if (i==0) {
      fprintf(fp,"set term png size 900,900\n");
      fprintf(fp,"set output \"%02d.png\"\n",idx);
    } else {
      fprintf(fp,"set term postscript eps enhanced color\n");
      fprintf(fp,"set output \"%02d.eps\"\n",idx);
    }
    fprintf(fp,"set size ratio -1\n");
    fprintf(fp,"set xrange [%d:%d]\n",roi[0],roi[1]);
    fprintf(fp,"set yrange [%d:%d]\n",roi[2],roi[3]);
    fprintf(fp,"set xlabel \"x [m]\"\n");
    fprintf(fp,"set ylabel \"z [m]\"\n");
    fprintf(fp,"plot \"%02d.txt\" using 1:2 lc rgb \"#FF0000\" title 'Ground Truth' w lines,",idx);
    fprintf(fp,"\"%02d.txt\" using 3:4 lc rgb \"#0000FF\" title 'Visual Odometry' w lines,",idx);
    fprintf(fp,"\"< head -1 %02d.txt\" using 1:2 lc rgb \"#000000\" pt 4 ps 1 lw 2 title 'Sequence Start' w points\n",idx);

    // close file
    fclose(fp);

    // run gnuplot => create png + eps
    sprintf(command,"cd %s; gnuplot %s",dir.c_str(),file_name);
    system(command);
  }

  // create pdf and crop
  sprintf(command,"cd %s; ps2pdf %02d.eps %02d_large.pdf",dir.c_str(),idx,idx);
  system(command);
  sprintf(command,"cd %s; pdfcrop %02d_large.pdf %02d.pdf",dir.c_str(),idx,idx);
  system(command);
  sprintf(command,"cd %s; rm %02d_large.pdf",dir.c_str(),idx);
  system(command);
}
