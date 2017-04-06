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
  //std::cout << "x_: " << x_ << std::endl;
  VectorXd y = z - H_ * x_;
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

  float ro = z[0];
  float phi = z[1];
  float nu = z[2];

  //MatrixXd Pj;
  //VectorXd x_temp;
  //Pj = MatrixXd(4, 3);
  //x_temp = VectorXd(3);
  //x_temp << ro, phi, nu;
  ////ekf_.x_ << ro * cos(phi), ro * sin(phi), 0, 0; // nu * cos(phi), nu * sin(phi);
  //Pj << cos(phi), -ro*sin(phi), 0,
  //  sin(phi), ro*cos(phi), 0,
  //  0, -nu*sin(phi), cos(phi),
  //  0, nu*cos(phi), sin(phi);
  //x_ = Pj * x_temp;

  x_ << ro * cos(phi), ro * sin(phi), nu * cos(phi), nu * sin(phi);

  H_ = tools.CalculateJacobian(x_);
  Update(z);

}
