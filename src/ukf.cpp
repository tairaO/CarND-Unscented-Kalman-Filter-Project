#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 6.0;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // time when the state is true, in us
  time_us_ = 0;

  // State dimension
  n_x_ = x_.size();

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  //set weights
  weights_.setConstant(weights_.size(), 1.0 / (2.0 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_+n_aug_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // the numbers that store measurements results
  n_meas_ = 0;
  n_over_ = 0;

  is_initialized_ = false;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if(!is_initialized_){
    // first measurement
    cout << "UKF: " << endl;
    // Set the initial mean and sigma
    x_ = VectorXd(n_x_);
    P_ = MatrixXd(n_x_, n_x_);
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // initialize state
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // initialize state from polar measurements.
      double rho_0 = meas_package.raw_measurements_[0];
      double phi_0 = meas_package.raw_measurements_[1];

      x_ << rho_0 * cos(phi_0), rho_0 * sin(phi_0), 0, 0, 0;
    }

    // set time
    time_us_ = meas_package.timestamp_;

    // finish initializing
    is_initialized_ = true;
    return;
  }
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // compute time step
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//[seconds]
 	time_us_ = meas_package.timestamp_;
  Prediction(dt);
 /*****************************************************************************
  *  Update
  ****************************************************************************/
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER)&& use_laser_) {
    UpdateLidar(meas_package);
  }
  else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)&&use_radar_) {
    UpdateRadar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  // print NIS
  cout << "measurements:" << n_meas_ << "  over 0.95 line:" << n_over_ << "  ratio(ideal:0.05):"<< n_over_/double(n_meas_) << endl;
  cout << "\n" << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  /*****************************************************
    *  Predict Sigma points
    *****************************************************/

  //create augmented mean vector and augumented state covariance
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augumented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state and augmented covariance matrix
  x_aug.fill(0); // Initialize
  x_aug.head(n_x_) = x_;

  P_aug.fill(0); // Initialize
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A_aug(n_aug_, n_aug_);
  A_aug = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for(int i=0; i < n_aug_; i++){
      Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }

  /*****************************************************
   *  Predict Sigma points
   *****************************************************/
   //extract values for better readability
   double p_x;
   double p_y;
   double v;
   double yaw;
   double yawd;
   double nu_a;
   double nu_yawdd;
   double dt2;
   dt2 = delta_t *delta_t;


   //predict sigma points
   for(int i=0; i < Xsig_aug.cols() ;i++){
      // extract values for easily reading code
       p_x = Xsig_aug(0, i);
       p_y = Xsig_aug(1, i);
       v   = Xsig_aug(2, i);
       yaw = Xsig_aug(3, i);
       yawd     = Xsig_aug(4, i);
       nu_a     = Xsig_aug(5, i);
       nu_yawdd = Xsig_aug(6, i);

       //avoid division by zero
       //write predicted sigma points into right column
       if(fabs(yawd) < 0.001){
         Xsig_pred_(0, i) = p_x + v*cos(yaw)*delta_t + dt2*cos(yaw)*nu_a/2.0;
         Xsig_pred_(1, i) = p_y + v*sin(yaw)*delta_t + dt2*sin(yaw)*nu_a/2.0;

       }else{
         Xsig_pred_(0, i) = p_x + v/yawd*( sin(yaw+yawd*delta_t)-sin(yaw))  + dt2*cos(yaw)*nu_a/2.0;
         Xsig_pred_(1, i) = p_y + v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw))  + dt2*sin(yaw)*nu_a/2.0;
       }
       Xsig_pred_(2, i) = v    + 0               + delta_t*nu_a;
       Xsig_pred_(3, i) = yaw  + yawd*delta_t    + dt2*nu_yawdd/2.0;
       Xsig_pred_(4, i) = yawd + 0               + delta_t*nu_yawdd;

   }


  /*****************************************************
   *  Predict Means and Covariances
   *****************************************************/

   // compute mean
   x_.fill(0.0); //initialize
   for(int i=0; i < Xsig_pred_.cols(); i++){
     //predict state mean
     x_ += weights_(i) * Xsig_pred_.col(i);
   }

   //compute covariance
   P_.fill(0.0); //initialize
   for(int i=0; i < Xsig_pred_.cols(); i++){
     //predict state covariance matrix
     VectorXd  x_diff  = Xsig_pred_.col(i) - x_;
     x_diff(3) = fmod(x_diff(3)+M_PI, 2.0*M_PI) - M_PI;
     P_ += weights_(i) * x_diff * x_diff.transpose() ;

   }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  /*****************************************************
   *  Predict Measurement
   *****************************************************/
  // recover measurements
  VectorXd z = meas_package.raw_measurements_;

  // number of sensor values
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // sensor noise matrix R
  MatrixXd R(n_z, n_z);
  R.fill(0);
  R(0, 0) = std_laspx_*std_laspx_;
  R(1, 1) = std_laspy_*std_laspy_;

  // create vector for storing sensor difference
  VectorXd z_diff(n_z);

  //extract values for better readability
  double p_x;
  double p_y;

  //transform sigma points into measurement space
  for(int i=0; i < Xsig_pred_.cols() ;i++){
   p_x  = Xsig_pred_(0, i);
   p_y  = Xsig_pred_(1, i);

   Zsig(0, i) = p_x;
   Zsig(1, i) = p_y;
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i=0; i < Xsig_pred_.cols() ;i++){
   z_pred += weights_(i) * Zsig.col(i);
  }
  //calculate innovation covariance matrix S
  S.fill(0.0);
  for(int i=0; i < Xsig_pred_.cols() ;i++){
   z_diff = Zsig.col(i) - z_pred;
   S += weights_(i) * z_diff * z_diff.transpose() ;
  }
  S += R;

  /*****************************************************
   *  Update state
   *****************************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // create matrix for Kalman gain K
  MatrixXd K(n_z, n_x_);

  // create vector for inovation vector
  VectorXd y = z - z_pred;

  // create vector for storing state difference
  VectorXd x_diff(n_x_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < Xsig_pred_.cols() ;i++){
    x_diff = Xsig_pred_.col(i) - x_;
    z_diff = Zsig.col(i) - z_pred;
    x_diff(3) = fmod(x_diff(3)+M_PI, 2.0*M_PI) - M_PI;
    Tc += weights_(i) * x_diff * z_diff.transpose();
}
//compute Kalman gain K
MatrixXd Si = S.inverse();
K = Tc * Si;

//update state mean and covariance matrix
x_ += K * y;
P_ -= K * S * K.transpose();

// print NIS
double nis;
nis = y.transpose() * Si * y;
cout << "NIS_L: " << nis << endl;

// update evaluation
n_meas_ += 1;
if(nis >5.991){
  n_over_ += 1;
  }

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  /*****************************************************
   *  Predict Measurement
   *****************************************************/

  // recover measurements
  VectorXd z = meas_package.raw_measurements_;

  // number of sensor values
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // sensor noise matrix R
  MatrixXd R(n_z, n_z);
  R.fill(0);
  R(0, 0) = std_radr_*std_radr_;
  R(1, 1) = std_radphi_*std_radphi_;
  R(2, 2) =std_radrd_*std_radrd_;

  // create vector for storing sensor difference
  VectorXd z_diff(n_z);

  //extract values for better readability
  double p_x;
  double p_y;
  double v;
  double yaw;


  //transform sigma points into measurement space
  for(int i=0; i < Xsig_pred_.cols() ;i++){
    p_x = Xsig_pred_(0, i);
    p_y = Xsig_pred_(1, i);
    v   = Xsig_pred_(2, i);
    yaw = Xsig_pred_(3, i);

    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/Zsig(0, i);
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i=0; i < Xsig_pred_.cols() ;i++){
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.fill(0.0);
  for(int i=0; i < Xsig_pred_.cols() ;i++){
    z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = fmod(z_diff(1)+M_PI, 2.0*M_PI) - M_PI;
    S += weights_(i) * z_diff * z_diff.transpose() ;
  }
  S += R;


  /*****************************************************
   *  Update state
   *****************************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // create matrix for Kalman gain K
  MatrixXd K;

  // create vector for inovation vector
  VectorXd y = z - z_pred;
  y(1) = fmod(y(1)+M_PI, 2.0*M_PI) - M_PI;

  // create vector for storing state difference
  VectorXd x_diff(n_x_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < Xsig_pred_.cols() ;i++){
    x_diff = Xsig_pred_.col(i) - x_;
    z_diff = Zsig.col(i) - z_pred;
    x_diff(3) = fmod(x_diff(3)+M_PI, 2.0*M_PI) - M_PI;
    z_diff(1) = fmod(z_diff(1)+M_PI, 2.0*M_PI) - M_PI;
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //compute Kalman gain K
  MatrixXd Si = S.inverse();
  K = Tc * Si;

  //update state mean and covariance matrix
  x_ += K * y;
  P_ -= K * S * K.transpose();

  // print NIS
  double nis;
  nis = y.transpose() * Si * y;
  cout << "NIS_R: " << nis << endl;

  // update evaluation
  n_meas_ += 1;
  if(nis >7.815){
    n_over_ += 1;
    }


}
