#include "ukf.h"
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
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  //Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // weights_ of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
     double weight = 0.5/(n_aug_+lambda_);
	 weights_(i) = weight;
   }

  // State covariance matrix
  P_ << 0.5, 0, 0, 0, 0,
		0, 0.5, 0, 0, 0,
    	    0, 0, 0.5, 0, 0,
    	    0, 0, 0, 0.5, 0,
    	    0, 0, 0, 0, 0.5;

  // Sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);


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
   if (!is_initialized_) {
	   InitializeUKF(meas_package);
	   return;
   }

   // Check the data source
   is_radar_ =  meas_package.sensor_type_ == MeasurementPackage::RADAR;

   // return if 1) radar data and ignore radar 2)laser data and ignore laser
   if((is_radar_ and !use_radar_) or (!is_radar_ and !use_laser_)) return;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//delta_t - expressed in seconds
   time_us_ = meas_package.timestamp_;
   Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
   UpdateSensor(meas_package);

   // false by default, will be updated in UpdateSensor() function
   is_radar_ = false;
}

/**
 * Initialize UKF with first measurement
 * @param meas_package The first measurement data of either radar or laser
 */
void UKF::InitializeUKF(MeasurementPackage meas_package) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
  	  float rho_measured = meas_package.raw_measurements_[0];
  	  float phi_measured = meas_package.raw_measurements_[1];
  	  float px = rho_measured * cos(phi_measured);
  	  float py = rho_measured * sin(phi_measured);
  	  x_ << px, py, 0, 0, 0;
      time_us_ = meas_package.timestamp_;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    	  x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      time_us_ = meas_package.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
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
   AugmentedSigmaPoints();
   SigmaPointPrediction(delta_t);
   PredictMeanAndCovariance();
}


/**
 * Updates the state and the state covariance matrix using either a laser
 * measurement or a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateSensor(MeasurementPackage meas_package) {

   if(is_radar_) n_z_ = 3;
   else n_z_ = 2;

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z_);
   //measurement covariance matrix S
   MatrixXd S = MatrixXd::Zero(n_z_,n_z_);
   //create matrix for sigma points in measurement space
   MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

   if(is_radar_) PredictRadarMeasurement(&z_pred, &S, &Zsig);
   else PredictLidarMeasurement(&z_pred, &S, &Zsig);

   VectorXd z  = meas_package.raw_measurements_;
   UpdateState(&z, &z_pred, &S, &Zsig);

}

/**
 * Augments the state and the state covariance matrix to create
 * augmented points
 * @param none
 */
void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

	P_aug = MatrixXd::Zero(7, 7);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6,6) = std_yawdd_ * std_yawdd_;

	Xsig_aug_.col(0) = x_aug;
	MatrixXd A = P_aug.llt().matrixL();
	for(int i = 0; i < n_aug_; i++) {
		Xsig_aug_.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

}

/**
 * Predicts single points
 * @param delta time
 */
void UKF::SigmaPointPrediction(double delta_t) {

    for(int i = 0; i < 2 * n_aug_ + 1; i++) {
      double px = Xsig_aug_(0, i); // x pos
      double py = Xsig_aug_(1, i); // y pos
      double v = Xsig_aug_(2, i); // velocity
      double psi = Xsig_aug_(3, i); // yaw
      double psi_dot = Xsig_aug_(4, i); // yawd
      double v_a = Xsig_aug_(5, i); // longitudinal acceleration noise
      double v_dd = Xsig_aug_(6, i); // yaw acceleration noise

      if(fabs(psi_dot) <= 0.001) {
    	  // divides by zero
    	    Xsig_pred_.col(i) << px + v*cos(psi)*delta_t + 0.5*delta_t*delta_t*cos(psi)*v_a,
              py + v*sin(psi)*delta_t + 0.5*delta_t*delta_t*sin(psi)*v_a,
              v + delta_t*v_a,
              psi + 0.5*delta_t*delta_t*v_dd,
              psi_dot + delta_t*v_dd;
      }
      else{
    	    Xsig_pred_.col(i) << px + v/psi_dot*(sin(psi+psi_dot*delta_t) - sin(psi)) + 0.5*delta_t*delta_t*cos(psi)*v_a,
              py + v/psi_dot*(cos(psi) - cos(psi+psi_dot*delta_t)) + 0.5*delta_t*delta_t*sin(psi)*v_a,
              v + delta_t*v_a,
              psi + psi_dot*delta_t + 0.5*delta_t*delta_t*v_dd,
              psi_dot + delta_t*v_dd;
      }
    }

}

/**
 * Predicts the state and the state covariance matrix
 * @param none
 */
void UKF::PredictMeanAndCovariance() {
	  x_.fill(0.0);
	  for(int i = 0; i < n_x_; i++) {
		  x_(i) = Xsig_pred_.row(i).dot(weights_);
	  }

	  P_.fill(0.0);
	  int yaw_pos = 3;
	  for(int i = 0; i < 2*n_aug_+1; i++) {
	     VectorXd diff = Xsig_pred_.col(i) - x_;
	     NormalizeAngle(&diff, yaw_pos);
	     P_ = P_ + weights_(i) * (diff * diff.transpose());
	  }
}

/**
 * Transforms sigma points into lidar measurements space and calculates mean and
 * covariance matrix
 * @param predicated mean, predicted Covariance, sigma points
 */
void UKF::PredictLidarMeasurement(VectorXd* z_pred,  MatrixXd* S, MatrixXd* Zsig) {

	 //transform sigma points into measurement space
	 for(int i = 0; i < 2 * n_aug_ + 1; i++) {
	   double px = Xsig_pred_(0, i); // x pos
	   double py = Xsig_pred_(1, i); // y pos
	   Zsig->col(i) << px, py;
	 }

	 //calculate mean predicted measurement
	 z_pred->fill(0.0);
	 for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	   //VectorXd tmp = VectorXd()
	   *z_pred = *z_pred + weights_(i) * Zsig->col(i);
	 }

	 //calculate measurement covariance matrix S
	 for(int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd diff = Zsig->col(i) - *z_pred;
		*S = *S + weights_(i) * (diff * diff.transpose());
	 }
	 MatrixXd R = MatrixXd(n_z_,n_z_);
	 R <<    std_laspx_*std_laspx_, 0,
			 0, std_laspy_*std_laspy_;

	 *S = *S + R;
}


/**
 * Transforms sigma points into radar measurements space and calculates mean and
 * covariance matrix
 * @param predicated mean, predicted Covariance, sigma points
 */
void UKF::PredictRadarMeasurement(VectorXd* z_pred,  MatrixXd* S, MatrixXd* Zsig) {

	 //transform sigma points into measurement space
	 for(int i = 0; i < 2 * n_aug_ + 1; i++) {
	   double px = Xsig_pred_(0, i); // x pos
	   double py = Xsig_pred_(1, i); // y pos
	   double v = Xsig_pred_(2, i); // velocity
	   double psi = Xsig_pred_(3, i); // yaw
	   double psi_dot = Xsig_pred_(4, i); // yawd

	   double rho, phi, rho_dot;
	   rho = sqrt(px*px + py*py);
	   if(px != 0) phi = atan2(py, px);
	   else phi = M_PI/2;
	   if(rho == 0) rho_dot = 0;
	   else rho_dot = (px*cos(psi)*v+py*sin(psi)*v)/rho;

	   Zsig->col(i) << rho, phi, rho_dot;
	 }

	 //calculate mean predicted measurement
	 z_pred->fill(0.0);
	 for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	   *z_pred = *z_pred + weights_(i) * Zsig->col(i);
	 }

	 //calculate measurement covariance matrix S
     int yaw_pos = 1;
	 for(int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd diff = Zsig->col(i) - *z_pred;
	    NormalizeAngle(&diff, yaw_pos);
		*S = *S + weights_(i) * (diff * diff.transpose());
	 }
	 MatrixXd R = MatrixXd(n_z_,n_z_);
	 R <<    std_radr_*std_radr_, 0, 0,
			 0, std_radphi_*std_radphi_, 0,
			 0, 0,std_radrd_*std_radrd_;

	 *S = *S + R;

}

/**
 * Update state mean and covariance matrix
 * @param measurements, predicated mean, predicted Covariance, sigma points
 */
void UKF::UpdateState(VectorXd* z, VectorXd* z_pred, MatrixXd* S,  MatrixXd* Zsig) {

	 //create matrix for cross correlation Tc
	  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);
	  //calculate cross correlation matrix
	  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
		  //residual
		  VectorXd z_diff = Zsig->col(i) - *z_pred;
		  if(is_radar_) NormalizeAngle(&z_diff, 1);

		  // state difference
		  VectorXd x_diff = Xsig_pred_.col(i) - x_;
		  NormalizeAngle(&x_diff, 3);

	      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	  }
	  //calculate Kalman gain K;
	  MatrixXd K = MatrixXd(n_x_, n_z_);
	  K = Tc * S->inverse();
	  //update state mean and covariance matrix
	  x_ = x_ + K * (*z - *z_pred);
	  P_ = P_ - K * *S * K.transpose();
}

/**
 * Normalizes angles
 * @param angles, yaw position in state vector
 */
void UKF::NormalizeAngle(VectorXd* diff, int yaw_pos) {
    while ((*diff)(yaw_pos)> M_PI) (*diff)(yaw_pos)-=2.*M_PI;
    while ((*diff)(yaw_pos)<-M_PI) (*diff)(yaw_pos)+=2.*M_PI;
}
