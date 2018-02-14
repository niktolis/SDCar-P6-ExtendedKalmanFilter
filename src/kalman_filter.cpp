#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in) {
  x_ = x_in;
  P_ = P_in;
}

MatrixXd KalmanFilter::GetStateCovMat()
{
  return P_;
}


void KalmanFilter::SetStateVec(VectorXd &x_in) {
  x_ = x_in;
}

VectorXd KalmanFilter::GetStateVec() {
  return x_;
}

void KalmanFilter::Predict(MatrixXd &F_in, MatrixXd &Q_in) {


  x_ = F_in * x_;

  MatrixXd Ft = F_in.transpose();
  P_ = F_in * P_ * Ft + Q_in;
}

void KalmanFilter::UpdateKF(const VectorXd &z, const MatrixXd &H_, const MatrixXd &R_) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  
  // call common update part
  Update(y, H_, R_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &H_, const MatrixXd &R_) {
  
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
  Update(y, H_, R_);
}

void KalmanFilter::Update(const VectorXd &y, const MatrixXd &H_, const MatrixXd &R_) {
  
  
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
