#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth){
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  // * the estimation vector size should not be zero
  // * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
            || estimations.size() == 0){
                
    cout << "CalculateRMSE() - Error - Invalid estimation or ground truth data!" << endl;
    return rmse;
  }
  
  // accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); i++){
      
      VectorXd residual = estimations[i] - ground_truth[i];
      
      // coefficient-wise multiplication
      residual = residual.array()*residual.array();
      
      rmse += residual;
  }
  
  // calculate the mean
  rmse = rmse /  estimations.size();
  
  // calculate the squared root
  rmse = rmse.array().sqrt();
  
  // return the result
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
  
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // pre-compute a set of terms to avoid repeated calculation
  float denom1 = px * px + py * py;
  float denom2 = sqrt(denom1);
  float denom3 = denom1 * denom2;
  
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
