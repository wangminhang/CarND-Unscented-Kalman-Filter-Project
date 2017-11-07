#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // false if class is first declared
  is_initialized_ = false;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.10;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.4;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = n_x_ + 2;

  //set augmented dimension
  n_sig_ = 2 * n_aug_ + 1;

  //create matrix for sigma points in measurement space of radar [r phi r_dot]
  Zsig_radar_ = MatrixXd::Zero(3, n_sig_);

  //create matrix for sigma points in measurement space of lidar [px py]
  Zsig_lidar_ = MatrixXd::Zero(2, n_sig_);

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // initial augmented state vector
  x_aug_ = VectorXd(7);

  // initial covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  // set initial augmented covariance matrix
  P_aug_ = MatrixXd::Zero(7, 7);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_,n_x_) = std_a_*std_a_;
  P_aug_(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //create sigma point matrix
  Xsig_aug_ = MatrixXd::Zero(n_aug_, n_sig_);

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_);

  // Initialize weights
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // Lidar measurement noise covariance matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

  // Radar measurement noise covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} m_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  /**
  *  Initialization
  */
  if (!is_initialized_)
  {
    double px = 0;
    double py = 0;
    double v = 0;


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
    }

    // If initial values are < 0.0001 px is set to default 0.1
    if(fabs(px) < 0.00001 && fabs(py) < 0.00001){
      px = 0.1;
      py = 0.1;
    }

    x_ << px, py, v, 0, 0;

    P_ << 0.15, 0, 0, 0, 0,
            0, 0.15, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    // timestamp
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
  *  Predict
  */
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  while (dt > 0.1)
  {
    const double delta_t = 0.05;
    Prediction(delta_t);
    dt -= delta_t;
  }
  Prediction(dt);

  /**
  *  Update
  */
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {


  /**
  * Generate Sigma Points
  */

  //create augmented mean state
  //augmented state vector made up of x_ (px,py,v,yaw,yawd) and 2 more points representing process noise
  //(acceleration noise, rate of change of yaw rate)
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //P_ is a 5x5 matrix for the original 5 states
  //P_aug is a 7x7 matrix to account for the 2 additional states added in x_aug_
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_,n_x_) = std_a_*std_a_;
  P_aug_(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  //15 sigma points are created by 2 x 7 states + 1
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++) //create other 14 sigma points
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  }

  /**
  * Predict Sigma Points
  */

  //n_sig_ is made up of 2 x (n_x + 2) + 1
  //2 times of 7 state variable + 1 to represent point at mean
  for (int i = 0; i< n_sig_; i++)
  {
    //extract values for better readability
    //run through column by column for all 15 sigma points
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else
    {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  /**
  * Predict Mean and Covariance
  */

  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++)
  {
    // state difference
    // taking individual sigma points - x_ (representing state mean)
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI)
        x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI)
        x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}



/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} m_pack
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  VectorXd z = meas_package.raw_measurements_;

  MatrixXd H_ = MatrixXd::Zero(2, 5);
  H_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;

  VectorXd z_pred = H_ * x_;
  VectorXd z_diff = z - z_pred;

  //update Kalman filter
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_lidar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * z_diff);
  P_ -= K * H_* P_;

  // Update NIS Lidar
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  VectorXd z = meas_package.raw_measurements_;

  /**
  * Use original predicted sigma points Xsig_pred_ and transfor them into measurement space r, phi, r_dot
  */

  for (int i = 0; i < n_sig_; i++) //2n+1 simga points
  {

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // convert to measurement model
    if(p_x == 0 && p_y == 0)
    {
      Zsig_radar_(0,i) = 0;   //r
      Zsig_radar_(1,i) = 0;   //phi
      Zsig_radar_(2,i) = 0;   //r_dot
    }
    else
    {
      Zsig_radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Zsig_radar_(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig_radar_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
  }

  /**
  *  predict mean and covariance for measurement model
  */

  VectorXd z_pred = VectorXd::Zero(3);
  for (int i=0; i < n_sig_; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig_radar_.col(i);
  }

  MatrixXd S = MatrixXd::Zero(3,3);
  for (int i = 0; i < n_sig_; i++) //2n+1 simga points
  {
    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI)
        z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI)
        z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance
  S = S + R_radar_;

  /**
  * UKF update
  1. Generate cross correlation matrix Tc
  2. Calculate Kalman gain taking Tc x S inverse
  */

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, 3); // 3 column for radar [r phi r_dot]

  //calculate cross correlation matrix
  for (int i = 0; i < n_sig_; i++)
  {
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred; // residual difference
    //angle normalization
    while (z_diff(1)> M_PI)
        z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI)
        z_diff(1) += 2.*M_PI;

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI)
        x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI)
        x_diff(3) += 2.*M_PI;

    //generate cross correlation matrix
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI)
      z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI)
      z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();

  // Update NIS Radar
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}