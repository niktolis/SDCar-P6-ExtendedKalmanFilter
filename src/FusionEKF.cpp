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
  
  cout << H_laser_.size() << endl;
  

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
              
              

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  
  /** 1. Initialize variables and matrices (x,F,H_laser, H_jacobian, P etc.)
   *  2. initialize the Kalman filter position vector with the first sensor measurements.
   *  3. modify the F and Q matrices prior to the prediction step based on the elapsed time between measurements.
   *  4. call the ipdate step for either the lidar or radar sensor measurement. 
   *     Because the update step for lidar and radar ae slightly different 
   *     there are different functions for updating lidar and radar.
  **/


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
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
        
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
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
      
      // set the state with the initial location and zero velocity
      x_ << px, py, vx, vy;
      
      
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      vx = 0;
      vy = 0;
      
    }
    
    // set the state with the initial location and zero velocity
    x_ << px, py, vx, vy;
    
    ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);    

    
    // save the initial timestamp to calculate dt for the next step.
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    
    
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;
   
   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;
   
   // Update the F Matrix with the new delta time
   ekf_.F_(0, 2) = dt;
   ekf_.F_(1, 3) = dt;
   
   ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
   
  cout << "Before predict!" << endl;
  ekf_.Predict();
  cout << "I Passed Predict" << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
    // Update the state matrix
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    // Update the covariance matrix
    ekf_.R_ = R_radar_;
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Update the state matrix
    ekf_.H_ = H_laser_;
    // Update the covariance matrix
    ekf_.R_ = R_laser_;
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
