#include "kalman_filter.h"
#include<iostream>

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
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Measurement residual.
  VectorXd y = z - H_ * x_;
  // Residual covariance.
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  // Optimal Kalman gain.
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // Posteriori state estimate.
  x_ = x_ + K * y;
  // Posteriori covariance estimate.
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Measurement residual.
  VectorXd y = z - H_ * x_;
  // Residual covariance.
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  // Optimal Kalman gain.
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // Posteriori state estimate.
  x_ = x_ + K*y;
  // Posteriori covariance estimate.
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
