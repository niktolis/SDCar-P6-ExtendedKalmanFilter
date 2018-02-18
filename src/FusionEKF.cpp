#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  F_ = MatrixXd(4, 4);
  P_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
    0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;

  F_ << 1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;

  P_ << 1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if(!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    VectorXd x_ = VectorXd(4);
    MatrixXd H_;

    double px, py, vx, vy;

    if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      // get the measurement data
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];


      px = rho * cos(phi);
      py = rho * sin(phi);
      vx = rho_dot * cos(phi);
      vy = rho_dot * sin(phi);

    } else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      vx = 0;
      vy = 0;
    }

    // set the state with the initial location and zero velocity
    x_ << px, py, vx, vy;

    ekf_.Init(x_, P_, F_, Q_, H_laser_);


    // save the initial timestamp to calculate dt for the next step.
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;



  } else {

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

     // new elapsed time. Time is measured in seconds.
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.Predict(dt, noise_ax_, noise_ay_);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      // Update the measurement jacobian matrix
      Hj_ = tools.CalculateJacobian(ekf_.GetStateVec());

      // Radar updates
      ekf_.Update(measurement_pack.raw_measurements_, R_radar_, Hj_);

    } else {

      // Laser updates
      ekf_.Update(measurement_pack.raw_measurements_, R_laser_);
    }

    // print the output
    cout << "x_ = " << ekf_.GetStateVec() << endl;
    cout << "P_ = " << ekf_.GetStateCovMat() << endl;
  }
}
