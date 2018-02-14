#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
private:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

public:
  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   */
  void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in);

  /**
   * Getter function for the state covariance
   * @return the state covariance matrix P_
   */
  Eigen::MatrixXd KalmanFilter::GetStateCovMat();

  /**
   * Setter function for the state vector
   * @param x_in State vector
   */
  void SetStateVec(Eigen::VectorXd & x_in);

  /**
  * Getter function for the state vector
  * @return the state vector x_
  */
  Eigen::VectorXd GetStateVec();


  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict(Eigen::MatrixXd &F_in, Eigen::MatrixXd &Q_in);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateKF(const Eigen::VectorXd &z, const Eigen::MatrixXd &H_, const Eigen::MatrixXd &R_);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z, const Eigen::MatrixXd &H_, const Eigen::MatrixXd &R_);

private:

  /**
   * Updates the state for KF and EKF given the vector y
   * @param y The result of the z_measured - z_predicted
   */
  void Update(const Eigen::VectorXd &y, const Eigen::MatrixXd &H_, const Eigen::MatrixXd &R_);

  /**
   * Normalizes the given angle in [-π, π)
   * @param phi The angle to be normalized in radians
   * @return the normalized angle in radians phi
   */
  double NormalizeAngle(double phi);

};

#endif /* KALMAN_FILTER_H_ */
