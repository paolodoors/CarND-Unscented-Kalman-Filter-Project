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
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8; // it was 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6; // it was 30

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  previous_timestamp_ = 0;
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;
  x_.fill(0.0);
  P_ = MatrixXd::Identity(n_x_, n_x_);
  //P_.fill(1.01);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ + 1);
  double w = lambda_ / (lambda_ + n_aug_);
  weights_(0) = w;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    w = 0.5 / (n_aug_ + lambda_);
    weights_(i) = w;
  }
  cout << "weights_:" << endl << weights_<< endl;

  NIS_radar_ = 0;
  NIS_laser_ = 0;
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
  cout << (meas_package.sensor_type_ == MeasurementPackage::RADAR ? "RADAR" : "LIDAR");
  cout << " measurement: ";
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_)
    || (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_)) {
    cout << "Discard!" << endl;
    return;
  } else {
    cout << endl << meas_package.raw_measurements_ << endl;
  }

  if (!is_initialized_) {
    // first measurement
    cout << "UKF:" << endl;

    previous_timestamp_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[2];
      float px = ro * cos(phi);
      float py = ro * sin(phi);
      //float vx = ro_dot * cos(phi);
      //float vy = ro_dot * sin(phi);

      // Discard bogus readings
      if (abs(ro) < 1e-10) {
        return;
      }
      x_(0) = px;
      x_(1) = py;
      is_initialized_ = true;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      /**
      Initialize state.
      */
      // Discard bogus readings
      if (abs(meas_package.raw_measurements_[0]) < 1e-10 || abs(meas_package.raw_measurements_[1]) < 1e-10) {
        return;
      }
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      is_initialized_ = true;
    }
    // print the output
    cout << "Init:" << endl;
    cout << "x_ = " << endl << x_ << endl;
    cout << "P_ = " << endl << P_ << endl;
    // done initializing, no need to predict or update
    return;
  }

  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    //Prediction(dt);
    //UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    Prediction(dt);
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
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
  cout << "Predict with delta_t: " << delta_t << endl;
  // augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;
  cout << "x_aug:" << endl << x_aug << endl;

  // augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.block(0, 0, n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = pow(std_a_, 2);
  P_aug((n_x_ + 1), (n_x_ + 1)) = pow(std_yawdd_, 2);
  cout << "P_aug:" << endl << P_aug << endl;

  // square root matrix P_aug
  MatrixXd A = P_aug.llt().matrixL();
  cout << "A:" << endl << A << endl;

  // X sigma augmented matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
  cout << "Xsig_aug:" << endl << Xsig_aug << endl;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid divizion by zero
    if (fabs(yawd) > 1e-3) {
      px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * pow(delta_t, 2) * cos(yaw);
    py_p = py_p + 0.5 * nu_a * pow(delta_t, 2) * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * pow(delta_t, 2);
    yawd_p = yawd_p + 0.5 * nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
  cout << "Xsig_pred_:" << endl << Xsig_pred_ << endl;

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  cout << ">>> x_: " << endl << x_ << endl;

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  cout << "Update LIDAR" << endl;
  // Set LASER dimension
  int n_z = 2;
  // sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }
  cout << "Zsig:" << endl << Zsig << endl;

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  cout << "z_pred:" << endl << z_pred << endl;

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  cout << "S:" << endl << S << endl;

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << pow(std_laspx_, 2), 0,
       0, pow(std_laspy_, 2);
  S = S + R;
  cout << "S:" << endl << S << endl;

  // matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  cout << "Tc:" << endl << Tc << endl;

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // compute NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
 cout << "End update" << endl;
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
}
