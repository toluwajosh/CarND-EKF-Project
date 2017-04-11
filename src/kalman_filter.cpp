#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

// initialize kalman filter
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
  /*
  predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
   /*
   update the state by using Kalman Filter equations
   */
  
  // prediction error
  VectorXd y = z - H_ * x_; 

  // perform rest of update
  UpdateOps(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /*
  update the state by using Extended Kalman Filter equations
  */

  // convert current state into polar coordinate
  float ro = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
  float phi = 0.0; // make phi 0.0 if too small
  if (fabs(x_[0]) > 0.001) {
    phi = atan2(x_[1], x_[0]);
  }
  float ro_dot = 0.0; // make ro_dot 0.0 if too small
  if (fabs(ro) > 0.001) {
    ro_dot = (x_[0] * x_[2] + x_[1] * x_[3]) / ro;
  }

  VectorXd hx(3);
  hx << ro, phi, ro_dot;
  
  // prediction error
  VectorXd y = z - hx;
  
  // perform rest of update
  UpdateOps(y);
}

void KalmanFilter::UpdateOps(const VectorXd &y) {
  /*
  perform rest of kalman filter update
  */

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si; // kalman gain

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
