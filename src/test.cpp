#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"
#include "tools.h"

using namespace std;


int main() {

  // Create a Kalman Filter instance
  UKF ukf;

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;
  MeasurementPackage meas_package;

  meas_package.

  //Call ProcessMeasurment(meas_package) for Kalman filter
  ukf.ProcessMeasurement(meas_package);

  return 0;
}
