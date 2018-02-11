#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  // test
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::UpdateKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  
  // call common update part
  Update(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  VectorXd z_pred = VectorXd(3);
  double px, py, vx, vy, rho, phi, rho_dot;
  
  px = x_(0);
  py = x_(1);
  vx = x_(2);
  vy = x_(3);
  
  rho = sqrt(px * px + py * py);
  phi = atan2(py, px);
  // check for division by zero
  if (fabs(rho) < 0.0001) {
	  std::cout << "UpdateEKF() - Error - Division by Zero!" << std::endl;
	  rho_dot = 0;
  }
  else
  {
	  rho_dot = (px * vx + py * vy) / rho;
  }
  
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  
  // Normalize the phi between [-π, π)
  y(1) = NormalizeAngle(y(1));
  
  // call common update part
  Update(y);
}

void KalmanFilter::Update(const VectorXd &y) {
  
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  // new estimate
  x_ = x_ + (K * y);
  int x_size = static_cast<int>(x_.size());
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


double KalmanFilter::NormalizeAngle(double phi){
    phi = fmod(phi + M_PI, 2 * M_PI);
    if (phi < 0)
        phi += 2 * M_PI;
    return phi - M_PI;
}
