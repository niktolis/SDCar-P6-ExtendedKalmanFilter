#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
private:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

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
   * @param P_in Initial state covariance matrix
   * @param F_in Initial state transition matrix
   * @param Q_in Initial process covariance matrix
   * @param H_in Initial measurement matrix
   */
  void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, 
            Eigen::MatrixXd &F_in, Eigen::MatrixXd &Q_in, 
            Eigen::MatrixXd &H_in);

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
   * @param dt Time between k and k+1 in s
   * @param noise_ax acceleration noise component in x axis
   * @param noise_ay accelaration noise component in y axis
   */
  void Predict(const double dt, const double noise_ax, const double noise_ay);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   * @param R_in The measurement covariance matrix
   */
  void Update(const Eigen::VectorXd &z, const Eigen::MatrixXd &R_in);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   * @param R_in The measurement covariance matrix
   * @param Hj_in The Jacobian measurement matrix
   */
  void Update(const Eigen::VectorXd &z, const Eigen::MatrixXd &R_in,
                 const Eigen::MatrixXd &Hj_in);

private:

  /**
   * The common part of the calculation of the KF and EKF
   * @param y The result of the z_measured - z_predicted
   * @param R_in The measurement covariance matrix
   * @param H_in The measurement matrix
   */
  void UpdateCommon(const Eigen::VectorXd &y, const Eigen::MatrixXd &R_in,
                    const Eigen::MatrixXd &H_in);

  /**
   * Normalizes the given angle in [-π, π)
   * @param phi The angle to be normalized in radians
   * @return the normalized angle in radians phi
   */
  double NormalizeAngle(double phi);

};

#endif /* KALMAN_FILTER_H_ */
