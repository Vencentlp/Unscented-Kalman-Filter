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
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
 
  is_initialized_ = false;
  
  //Initialize process covariance P_
  P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 1, 0,
	  0, 0, 0, 0, 1;

  // Initialize some parameters
  time_us_ = 0;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2 * n_aug_ + 1);
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (!is_initialized_) 
	{
		if (meas_package.sensor_type_ == meas_package.LASER)
		{
			double px = meas_package.raw_measurements_(0);
			double py = meas_package.raw_measurements_(1);
			x_ << px, py, 0, 0, 0;
		}
		else
		{
			double r = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			//double r_dot = meas_package.raw_measurements_(2);


			double px = r * cos(phi);
			double py = r * sin(phi);
			//double v = r_dot;
			x_ << px, py, 0, 0, 0;
		}
		is_initialized_ = true;	
		time_us_ = meas_package.timestamp_;

		//MatrixXd R_laser = MatrixXd(2, 2);
		//R_laser << std_laspx_ * std_laspx_, 0,
			//0, std_laspy_*std_laspy_;		
	}
	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	std::cout << "delta_t" << delta_t << std::endl;
	Prediction(delta_t);
	if (meas_package.sensor_type_ == meas_package.LASER)
	{
		UpdateLidar(meas_package);
	}
	if (meas_package.sensor_type_ == meas_package.RADAR)
	{
		UpdateRadar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	// Creat augmented mean state
	VectorXd x_aug = VectorXd(n_aug_);

	// Creat augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//Creat sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	
	//Creat augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//Creat state noise covariance matrix
	MatrixXd Q_ = MatrixXd(2, 2);
	Q_ << std_a_ * std_a_, 0,
		0, std_yawdd_*std_yawdd_;

	// Creat state covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug.bottomRightCorner(2, 2) = Q_;
	
	std::cout << "Paug" << P_aug << std::endl;
	//Creat square root matrix
	MatrixXd L = P_aug.llt().matrixL();
	std::cout << "L" << L << std::endl;
	
	//Creat sigma points
	Xsig_aug.col(0) = x_aug;
	
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_)*L.col(i);
		Xsig_aug.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_)*L.col(i);
	}
	
	
	//predict sigma points
	Xsig_pred_ = MatrixXd(5, 2 * n_aug_ + 1);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double px = Xsig_aug(0, i);
		double py = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);
		double px_p, py_p;
		if (abs(yawd) < 0.0001)
		{
			px_p = px + v * delta_t*cos(yaw);
			py_p = py + v * delta_t*sin(yaw);
		}
		else
		{
			px_p = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = py + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		}
		
		double v_p = v + delta_t * nu_a;
		double yaw_p = yaw + yawd * delta_t + 0.5*delta_t*delta_t*nu_yawdd;
		double yawd_p = yawd + delta_t * nu_yawdd;
		
		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
		//Xsig_pred_.col(i) << px_p, py_p, v_p, yaw_p, yawd_p;
	}
	
	
	// predict mean and covariance for sigma points

	weights_.fill(0.5 / (lambda_ + n_aug_));
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	
	x_.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	
	VectorXd x_diff = VectorXd(n_x_);
	MatrixXd P_diff = MatrixXd(n_x_, n_x_);

	P_.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_diff = Xsig_pred_.col(i) - x_;
		while (x_diff(3) < -M_PI) { x_diff(3) += 2 * M_PI; }
		while (x_diff(3) > M_PI) { x_diff(3) -= 2 * M_PI; }
		P_diff = weights_(i)*x_diff*x_diff.transpose();
		P_ = P_ + P_diff;
	}
	std::cout << "Predictx " << x_ << std::endl;
	std::cout << "PredictP " << P_ << std::endl;
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
	MatrixXd R_laser = MatrixXd(2, 2);
	R_laser << std_laspx_ * std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
	Zsig.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);

		Zsig(0, i) = px;
		Zsig(1, i) = py;
	}

	VectorXd z_pred = VectorXd(2);
	z_pred.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred += weights_(i) * Zsig.col(i);
	}

	MatrixXd S = MatrixXd(2, 2);
	S.fill(0);
	VectorXd z_diff = VectorXd(2);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_diff = Zsig.col(i) - z_pred;
		S += weights_(i)*z_diff*z_diff.transpose();
	}
	S += R_laser;

	//Update
	MatrixXd T = MatrixXd(5, 2);
	T.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;		
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		T += weights_(i)*x_diff*z_diff.transpose();
	}
	MatrixXd K = T * S.inverse();
	VectorXd z = VectorXd(2);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
	x_ = x_ + K * (z - z_pred);
	P_ = P_ - K * S*K.transpose();
	std::cout << "LaserT " << T << std::endl;
	std::cout << "LaserK " << K << std::endl;
	std::cout << "LaserP " << P_ << std::endl;
	std::cout << "Laserx" << x_ << endl;
	double NIS = (z - z_pred).transpose()*(S.inverse())*(z - z_pred);
	std::cout << "Laser NIS:" << NIS << std::endl;
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
	MatrixXd R_Radar = MatrixXd(3, 3);
	R_Radar << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	
	MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
	Zsig.fill(0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);
		double v1 = v * cos(yaw);
		double v2 = v * sin(yaw);

		Zsig(0, i) = sqrt(px*px + py * py);
		if (abs(px) > 0.0001)
		{
			Zsig(1, i) = atan2(py, px);
		}
		else
		{
			std::cout << "0-atan2" << Zsig(1, i) << std::endl;
		}
		if (Zsig(0, i) > 0.0001)
		{			
			Zsig(2, i) = (px*v1 + py * v2) / sqrt(px*px + py * py);
		}
		else
		{
			std::cout << "0-r_dot" << Zsig(2, i) << std::endl;			
		}
		
		
	}
	
	VectorXd z_pred = VectorXd(3);
	z_pred.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred += weights_(i) * Zsig.col(i);
		
	}
	while (z_pred(1) < -M_PI) { z_pred(1) += 2 * M_PI; }
	while (z_pred(1) > M_PI) { z_pred(1) -= 2 * M_PI; }
	
	MatrixXd S = MatrixXd(3, 3);
	S.fill(0);
	VectorXd z_diff = VectorXd(3);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_diff = Zsig.col(i) - z_pred;
		while (z_diff(1) < -M_PI) { z_diff(1) += 2 * M_PI; }
		while (z_diff(1) > M_PI) { z_diff(1) -= 2 * M_PI; }
		S += weights_(i)*z_diff*z_diff.transpose();
	}
	
	S += R_Radar;


	//Update
	MatrixXd T = MatrixXd(5, 3);
	T.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		T += weights_(i)*x_diff*z_diff.transpose();
	}
	
	MatrixXd K = T * S.inverse();
	VectorXd z = VectorXd(3);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
	z_diff = z - z_pred;
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S*K.transpose();
	std::cout << "RadarT " << T << std::endl;
	std::cout << "RadarK " << K << std::endl;
	std::cout << "RadarP " << P_ << std::endl;
	std::cout << "Radarx" << x_ << endl;
	double NIS = (z - z_pred).transpose()*(S.inverse())*(z - z_pred);
	std::cout << "Radar NIS:" << NIS << std::endl;
}
