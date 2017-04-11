#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// Constructor.
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initialize matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // initialize jacobian matrix - Radar. Will be later updated
  Hj_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  // process noise
  noise_ax = 9;
  noise_ay = 9;
}

// Destructor.
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {

    // Initialize the state ekf_.x_
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // prediction covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.Q_ = MatrixXd(4, 4);

    // initialize transition matrix
    ekf_.F_ << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

    // initialize prediction covariance matrix
    ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      // convert radar measurements from polar to cartesian
      // Set state ekf_.x_ to the first measurement.
      double ro = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double ro_dot = measurement_pack.raw_measurements_[2];

      ekf_.x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Set state ekf_.x_ to the first measurement.
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0.0, 0.0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // update F matrix
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Update the process covariance matrix Q
  // noise_ax = 9 and noise_ay = 9
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
    0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
    dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
    0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

  // predict
  ekf_.Predict();

  /*****************************************************************************
   *  Update: the state and covariance matrices.
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
    // Radar updates
    Tools tools;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    
  } else {

    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
