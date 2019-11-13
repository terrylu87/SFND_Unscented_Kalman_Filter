#include "ukf.h"
#include "Eigen/Dense"
#include <vector>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const float EPSILON = 0.0001;
//#define USE_UKF_FOR_LIDAR


/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
    P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

    n_x_ = 5;

    n_aug_ = 7;

    lambda_ = 3 - n_x_;

    weights_ =  VectorXd(2*n_aug_+1);
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }

    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;

    R_lidar_ = MatrixXd(2, 2);
    R_lidar_ << std_laspx_*std_laspx_,0,
        0,std_laspy_*std_laspy_;

    H_lidar_ = MatrixXd(2,5);
    H_lidar_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

    Q_ = MatrixXd(2,2);
    Q_ << std_a_*std_a_,0,
        0,std_yawdd_*std_yawdd_;
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
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR
        && use_radar_ == false) {
        return;
    }
    if (meas_package.sensor_type_ == MeasurementPackage::LASER
        && use_laser_ == false) {
        return;
    }
    if (!is_initialized_) {
        /**
           TODO:
           * Initialize the state x_ with the first measurement.
           * Create the covariance matrix.
           * Remember: you'll need to convert radar from polar to cartesian coordinates.
           */
        // first measurement
        //cout << "UKF: " << endl;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
               Convert radar from polar to cartesian coordinates and initialize state.
            */
            float px = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
            float py = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
            float v = sqrt(px*px+py*py);
            x_ << px, py, v, 0, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
               Initialize state.
            */
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }

        if (fabs(x_(0)) < EPSILON and fabs(x_(1)) < EPSILON){
            x_(0) = EPSILON;
            x_(1) = EPSILON;
        }

        previous_timestamp_ = meas_package.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }


    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(dt);


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        UpdateRadar(meas_package);
    } else {
        // Laser updates
        UpdateLidar(meas_package);
    }

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
    // 1. generate sigma points
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(7, 7);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug.topLeftCorner(5,5) = P_;
    P_aug.block(5,5,2,2) = Q_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    // 2. predict sigma points
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        const double p_x = Xsig_aug(0,i);
        const double p_y = Xsig_aug(1,i);
        const double v = Xsig_aug(2,i);
        const double yaw = Xsig_aug(3,i);
        const double yawd = Xsig_aug(4,i);
        const double nu_a = Xsig_aug(5,i);
        const double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
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
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
    Xsig_pred_ = Xsig_pred;
    // 3. predict mean and covariance
    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);
    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);
    //predicted state mean
    x.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x = x+ weights_(i) * Xsig_pred.col(i);
    }

    //predicted state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P = P + weights_(i) * x_diff * x_diff.transpose();
    }
    P_ = P;
    x_ = x;

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

    // use common Kalman Filter

    VectorXd z = meas_package.raw_measurements_;
#ifdef USE_UKF_FOR_LIDAR
    // 1. Predit measurement
    int n_z = 2;
    MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    S = S + R_lidar_;

    // 2. Update state
    // Incoming radar measurement
    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
#else
	VectorXd z_pred = H_lidar_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_lidar_.transpose();
	MatrixXd S = H_lidar_ * P_ * Ht + R_lidar_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_lidar_) * P_;
#endif

    ////residual
    VectorXd z_diff = z - z_pred;

#ifdef USE_UKF_FOR_LIDAR
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
#endif

    //NIS Lidar Update
    NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
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
    
    int n_z = 3;
    int n_aug = n_aug_;
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    // 1. predict measurement
    MatrixXd Xsig_pred = Xsig_pred_;
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v  = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    S = S + R_radar_;

    // 2. update state
    VectorXd z = meas_package.raw_measurements_;

    //calculate cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    // print the output
    //cout << "x_ = " << x_ << endl;
    //cout << "P_ = " << P_ << endl;

    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
