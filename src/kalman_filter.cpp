#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

double KalmanFilter::normalize(double x){
  if (x < 0)
    x += 360;
  x = fmod(x + 180,360);
  return x - 180;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  float pxpysqr = px*px + py*py;
  float pxpysqrt = sqrt(pxpysqr);
  float rho = pxpysqrt;
  float phi = 0 ;
  float rho_dot = 0;
  if(vx !=0){
    phi = atan(vy/vx);
  }
  if(pxpysqrt > 0.01f){
    rho_dot = (px*vx + py*vy)/pxpysqrt;
  }
  VectorXd z_pred = VectorXd(3);

  z_pred << rho,phi,rho_dot;

  VectorXd y = z - z_pred;

  y(1) = normalize(y(1));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;


  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}

//int main(int argc, char* argv[]){
//  KalmanFilter kf = KalmanFilter();
//  cout<< kf.normalize(0) <<endl;
//  cout<< kf.normalize(90)<<endl;
//  cout<< kf.normalize(180)<<endl;
//  cout<< kf.normalize(270)<<endl;
//  cout<< kf.normalize(360)<<endl;
//  return 0;
//}
