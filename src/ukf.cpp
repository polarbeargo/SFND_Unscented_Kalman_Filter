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

  // Process noi5e standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  // set augmented dimension
  n_aug_ = 7;

  // expressed in a variable improve performance
  n_sig_ = 2 * n_aug_ + 1;

  // Set initialize to false as default
  is_initialized_ = false;

  // Time when the state is true, in us
  time_us_ = 0;
  n_x_ = 5;

  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(5, n_sig_);

  // set weights
  weights_ = VectorXd(n_sig_);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  // Generating Sigma Points
  if (!is_initialized_)
  {

    switch (meas_package.sensor_type_)
    {
    case MeasurementPackage::LASER:
      x_(0) = meas_package.raw_measurements_(0); // p_x
      x_(1) = meas_package.raw_measurements_(1); // p_y
      break;

    case MeasurementPackage::RADAR:
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);
      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = rho_d;
      break;
    }

    P_.setIdentity();
    P_(2, 2) = 100.;
    P_(3, 3) = 10.;

    is_initialized_ = true;
    return;
  }
  // Predict
  double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  switch (meas_package.sensor_type_)
  {
  case MeasurementPackage::LASER:
    UpdateLidar(meas_package);
    break;
  case MeasurementPackage::RADAR:
    UpdateRadar(meas_package);
    break;
  }
}

void UKF::Prediction(double delta_t)
{
  // /**
  //  * TODO: Complete this function! Estimate the object's location.
  //  * Modify the state vector, x_. Predict sigma points, the state,
  //  * and the state covariance matrix.
  //  */
  // Augmentation
  VectorXd x_aug = VectorXd(7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd P_aug = MatrixXd(7, 7);

  x_aug.head(5) = x_;
  x_aug[5] = 0;
  x_aug[6] = 0;

  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // Create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  MatrixXd L = P_aug.llt().matrixL();
  Xsig_aug.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    Xsig_aug.col(i) = x_aug;
  }

  Xsig_aug.block<7, 7>(0, 1) += sqrt(lambda_ + n_aug_) * L;
  Xsig_aug.block<7, 7>(0, n_aug_ + 1) -= sqrt(lambda_ + n_aug_) * L;

  // Sigma Point Prediction

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(n_x_, i);
    double nu_yawdd = Xsig_aug(n_x_ + 1, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
  // Predict Mean and Covariance
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
  {
    weights_(i) = 1 / (2 * lambda_ + 2 * n_aug_);
  } // set weights
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    x_ += weights_(i) * Xsig_pred_.col(i);
  } // predict state mean
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff, 3);
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  VectorXd z = meas_package.raw_measurements_;
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  MatrixXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Zsig.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }
  z_pred.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  } // mean predicted measurement

  S.setZero();
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff, 3);

    Tc += weights_(i) * x_diff * z_diff.transpose();
  } // innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
  S += R; // add measurement noise covariance matrix

  MatrixXd K = Tc * S.inverse(); // Kalman gain K;
  VectorXd z_diff = z - z_pred;  // residual

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // Predict Radar Measurement
  VectorXd z = meas_package.raw_measurements_;
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  MatrixXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Zsig.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);
  }
  z_pred.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  } // mean predicted measurement

  S.setZero();
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff, 1);
    S += weights_(i) * z_diff * z_diff.transpose();
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff, 3);

    Tc += weights_(i) * x_diff * z_diff.transpose();
  } // innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  S += R; // add measurement noise covariance matrix

  MatrixXd K = Tc * S.inverse(); // Kalman gain K;
  VectorXd z_diff = z - z_pred;  // residual
  NormalizeAngle(z_diff, 1);
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}
/**
 *   normalize angles between [-M_PI, M_PI] interval.
 */
void UKF::NormalizeAngle(Eigen::VectorXd &vector, int index)
{
  while (vector(index) > M_PI)
    vector(index) -= 2. * M_PI;
  while (vector(index) < -M_PI)
    vector(index) += 2. * M_PI;
}