#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  std::cout << "Last Ground truth:" << ground_truth[1223] << endl;
  std::cout << "Last Estimation:" << estimations[1223] << endl;
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;


  // Sum squared difference between estimation and ground truth.
  for (int i=0; i < ground_truth.size(); i++){
    // Calculate error.
    VectorXd error = ground_truth[i] - estimations[i];

    // Square the error.
    error = error.array() * error.array();

    // Sum the error.
    rmse += error;
  }

  // Find mean error.
  rmse = rmse/ground_truth.size();

  // Find square root of mean squared error.
  rmse = rmse.array().sqrt();

  return rmse;
  /**
  TODO:
    * Calculate the RMSE here.
  */
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  // Set state parameters to help readability.
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];

  // Compute common sets used in the jacobian function.
  float xysq = px*px + py*py;
  float xysqrt = sqrt(xysq);
  float xy32 = xysq * xysqrt;

  MatrixXd Hj(3,4);

  Hj << px/xysqrt, py/xysqrt, 0, 0,
        -py/xysq, px/xysq, 0, 0,
        py*(vx*py - vy*px)/xy32, px*(px*vy - py*vx)/xy32, px/xysqrt, py/xysqrt;
  return Hj;

  
  /**
  TODO:
    * Calculate a Jacobian here.s
  */
}
