#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {

  product_sum_ = VectorXd(4);
  product_sum_ << 0.0, 0.0, 0.0, 0.0;

}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {


  VectorXd rmse;

  // check the validity of the following inputs:
  // * the estimation vector size should not be zero
  // * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0) {

    cout << "CalculateRMSE() - Error - Invalid estimation or ground truth data!"
      << endl;
  } else {

    unsigned int n = static_cast<unsigned int>(estimations.size());

    VectorXd residual = estimations[n - 1] - ground_truth[n - 1];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();

    // accumulate squared residuals
    product_sum_ += residual;

  }

  // calculate the mean
  rmse = product_sum_ / static_cast<int>(estimations.size());

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3, 4);

  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  double d1 = px * px + py * py;
  double d2 = sqrt(d1);
  double d3 = d1 * d2;

  // check for division by zero
  if(fabs(d1) < 0.0001) {
    cout << "CalculateJacobian() - Error - Division by Zero!" << endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px / d2), (py / d2), 0, 0,
    -(py / d1), (px / d1), 0, 0,
    py * (py * vx - px * vy) / d3, px * (px * vy - py * vx) / d3, px / d2, py / d2;

  return Hj;

}
