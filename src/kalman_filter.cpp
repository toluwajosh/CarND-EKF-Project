#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  std::cout << "X predicted: " <<x_<< std::endl;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  //std::cout << "z: " << z << std::endl;
  //std::cout << "H_: " << H_ << std::endl;
  std::cout << "F matrix: " << F_ << std::endl;
  VectorXd y = z - H_ * x_;
  std::cout << "y: " << y << std::endl;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  //std::cout << x_ << std::endl;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  Tools tools;

  /////////////////////////////////////////////////////////////////
  //float ro;
  //float phi;
  //float nu;

  //ro = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
  //phi = atan2(x_[1], x_[0]);
  //nu = (x_[0] * x_[2] + x_[1] * x_[3]) / ro;

  //VectorXd xx;
  //xx = VectorXd(3);
  //xx << ro, phi, nu;
  /////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  //float ro=z[0];
  //float phi=z[1];
  //float nu=z[2];
  //x_ << ro * cos(phi), ro * sin(phi), nu * cos(phi), nu * sin(phi);
  //H_ = tools.CalculateJacobian(x_);
  //////////////////////////////////////////////////////////////////

  // the update equations ////////////////////////////////////////////////////
  VectorXd hx = H_ * x_;
  VectorXd y = z - hx;
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
