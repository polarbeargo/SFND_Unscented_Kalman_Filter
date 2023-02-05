#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  // Set initialize to false as default
  is_initialized_ = false;

  // Time when the state is true, in us
  time_us_ = 0.0;

  // Define spreading parameter
  lambda_ = 3 - n_x_;

  // set augmented dimension
  n_aug_ = 7;

  // set measurement dimension
  int n_z = 2;

  // expressed in a variable improve performance
  n_sig_ = 2 * n_aug_ + 1;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // create vector for weights
  weights_ = VectorXd(n_sig_);

  // Initialize measurement noice covarieance matrix
  R_radar_ = MatrixXd(n_z + 1, n_z + 1);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(n_z, n_z);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
    time_us_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {

      /**
      Convert radar from polar to cartesian coordinates
      */

      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      float px = rho * cos(phi);
      float py = rho * sin(phi);

      x_ << px, py, 0, 0, 0;
    }
    else
    {
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];

      x_ << px, py, 0, 0, 0;
    }

    // Saving First timestamp
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  // switch between lidar and radar measurements
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t)
{
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}