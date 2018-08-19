#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse_(4);
  rmse_ << 0, 0, 0, 0;
  unsigned int est_num_ =  estimations.size();

    if(est_num_ != ground_truth.size() || est_num_==0) {
      cout << "Invalid estimation or ground_truth data" << endl;
      return rmse_;
    }

    // accumulate residuals
    for(unsigned int i=0; i < est_num_; ++i) {
      // calculate squared error
      VectorXd residual_ = estimations[i] - ground_truth[i];
      residual_ = residual_.array()*residual_.array();
      rmse_ += residual_;
    }

    // calculate the mean
    rmse_ = rmse_/est_num_;
    // calculate the root
    rmse_ = rmse_.array().sqrt();

    return rmse_;
}
