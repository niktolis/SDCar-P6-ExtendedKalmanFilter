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
                
    cout << "CalculateRMSE() - Error - Invalid estimation or ground truth data!" << endl;
  }
  else {

    unsigned int n = estimations.size();
    cout << "px_est = " << estimations[n - 1](0) << "\t" << "px_gt = " << ground_truth[n - 1](0) << endl;
    cout << "py_est = " << estimations[n - 1](1) << "\t" << "py_gt = " << ground_truth[n - 1](1) << endl;
    cout << "vx_est = " << estimations[n - 1](2) << "\t" << "vx_gt = " << ground_truth[n - 1](2) << endl;
    cout << "vy_est = " << estimations[n - 1](3) << "\t" << "vy_gt = " << ground_truth[n - 1](3) << endl;

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

  MatrixXd Hj(3,4);
  
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  // pre-compute a set of terms to avoid repeated calculation
  double denom1 = px * px + py * py;
  double denom2 = sqrt(denom1);
  double denom3 = denom1 * denom2;
  
  // check for division by zero
  if(fabs(denom1) < 0.0001){
      cout << "CalculateJacobian() - Error - Division by Zero!" << endl;
      return Hj;
  }
  
  // compute the Jacobian matrix
  Hj << (px/denom2), (py/denom2), 0, 0,
        -(py/denom1), (px/denom1), 0, 0,
        py*(py*vx - px*vy)/denom3, px*(px*vy - py*vx)/denom3, px/denom2, py/denom2;
  
  return Hj;
  
}
