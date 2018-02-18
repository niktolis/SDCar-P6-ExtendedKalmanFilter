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
                        MatrixXd &Q_in, MatrixXd &H_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  Q_ = Q_in;
  H_ = H_in;
}

MatrixXd KalmanFilter::GetStateCovMat() {
  return P_;
}


void KalmanFilter::SetStateVec(VectorXd &x_in) {
  x_ = x_in;
}

VectorXd KalmanFilter::GetStateVec() {
  return x_;
}

void KalmanFilter::Predict(const double dt, const double noise_ax,
                           const double noise_ay) {

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  // Update the values of the state transition matrix based on
  // the new elapset time
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  // Update the process covariance matrix
  Q_ << (dt_4*noise_ax) / 4, 0, (dt_3*noise_ax) / 2, 0,
    0, (dt_4*noise_ay) / 4, 0, (dt_3*noise_ay) / 2,
    (dt_3*noise_ax) / 2, 0, dt_2*noise_ax, 0,
    0, (dt_3*noise_ay) / 2, 0, dt_2*noise_ay;

  x_ = F_ * x_;

  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &R_in) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  // call common update part
  UpdateCommon(y, R_in, H_);
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &R_in, 
                          const MatrixXd &Hj_in) {

  VectorXd z_pred = VectorXd(3);
  double px, py, vx, vy, rho, phi, rho_dot;

  px = x_(0);
  py = x_(1);
  vx = x_(2);
  vy = x_(3);

  rho = sqrt(px * px + py * py);
  phi = atan2(py, px);
  // check for division by zero
  if(fabs(rho) < 0.0001) {
    std::cout << "UpdateEKF() - Error - Division by Zero!" << std::endl;
    rho_dot = 0;
  } else {
    rho_dot = (px * vx + py * vy) / rho;
  }

  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;

  // Normalize the phi between [-π, π)
  y(1) = NormalizeAngle(y(1));

  // call common update part
  UpdateCommon(y, R_in, Hj_in);
}

void KalmanFilter::UpdateCommon(const VectorXd &y, const MatrixXd &R_,
                                const MatrixXd &H_in) {


  MatrixXd Ht = H_in.transpose();
  MatrixXd S = H_in * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  // new estimate
  x_ = x_ + (K * y);
  int x_size = static_cast<int>(x_.size());
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_in) * P_;
}


double KalmanFilter::NormalizeAngle(double phi) {
  phi = fmod(phi + M_PI, 2 * M_PI);
  if(phi < 0)
    phi += 2 * M_PI;
  return phi - M_PI;
}
