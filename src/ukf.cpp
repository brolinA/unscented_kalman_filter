#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.0;
  
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

  is_initialized_ = false;
  
  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, (2*n_aug_)+1);

  // add measurement noise covariance matrix
  noise_radar_ = MatrixXd(3, 3);
  noise_radar_.fill(0.0);
  noise_radar_(0, 0) = std_radr_ * std_radr_;
  noise_radar_(1, 1) = std_radphi_ * std_radphi_;
  noise_radar_(2, 2) = std_radrd_ * std_radrd_;


  noise_lidar_ = MatrixXd(2, 2);
  noise_lidar_.fill(0.0);
  noise_lidar_(0, 0) = std_laspx_ * std_laspx_;
  noise_lidar_(1, 1) = std_laspy_ * std_laspy_;

  //filling the weight vector with the values.
  weights_ = VectorXd((2*n_aug_)+1);
  weights_.fill(0.5/(lambda_+n_aug_));
  weights_(0) = lambda_/(lambda_+n_aug_);

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            0.1,
            0.1,
            0.1;

      double x_variance = std_laspx_*std_laspx_;
      double y_variance = std_laspy_*std_laspy_;
      
      P_ << x_variance,     0,        0, 0, 0,
                0,      y_variance,   0, 0, 0,
                0,          0,        1, 0, 0,
                0,          0,        0, 1, 0,
                0,          0,        0, 0, 1;

      std::cout << "Initialized with LIDAR data" << std::endl;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double range = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      double theta_dot = meas_package.raw_measurements_(2);

      x_ << range*cos(theta),
            range*sin(theta),
            0,
            0,
            0;

      double range_variance = std_radr_*std_radr_;
      double theta_variance = std_radrd_*std_radrd_;

      P_ << range_variance,       0,                0,           0,          0,
                 0,         range_variance,         0,           0,          0,
                 0,               0,        theta_variance,      0,          0,
                 0,               0,                0,      std_radphi_,     0,
                 0,               0,                0,           0,      std_radphi_;

      std::cout << "Initialized with RADAR data" << std::endl;
    }
    else
    {
      std::cout << "This UKF is no designed for data from " << meas_package.sensor_type_ << " sensor type. Input a valid sensor data." << std::endl;
      return;
    }

    // P_ = MatrixXd::Identity(5,5);

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    std::cout << "Filter initliazed successfully" << std::endl;
    return;
  }
  
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    UpdateLidar(meas_package);
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);

}

void UKF::Prediction(double delta_t) 
{

  //1st step is to create sigma points
  
  //creating augmented state vector x_aug
  Eigen::VectorXd x_aug = VectorXd(n_aug_);
  //x_aug << x_, 0.0, 0.0;
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  //creating augmented covariance matrix
  Eigen::MatrixXd p_aug = MatrixXd(n_aug_, n_aug_);
  p_aug.fill(0.0);
  p_aug.topLeftCorner(n_x_,n_x_) = P_;
  p_aug(n_x_, n_x_) = std_a_*std_a_;
  p_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

  //root of the covariance matrix
  Eigen::MatrixXd root_p = p_aug.llt().matrixL();

  //create the sigma points matrix.
  Eigen::MatrixXd x_sig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  x_sig_aug.fill(0.0);
  x_sig_aug.col(0) = x_aug;

  for (int i = 0; i <n_aug_; i++)
  {
    x_sig_aug.col(i+1) = x_aug + (sqrt(lambda_+n_aug_)*root_p.col(i));
    x_sig_aug.col(i+1+n_aug_) = x_aug - (sqrt(lambda_+n_aug_)*root_p.col(i));
  }

  //predicting sigma points
  for (int i = 0; i <= 2*n_aug_; i++)
  {
    double px = x_sig_aug(0,i);
    double py = x_sig_aug(1,i);
    double v = x_sig_aug(2,i);
    double psy = x_sig_aug(3,i);
    double psy_dot = x_sig_aug(4,i);
    double nu_a = x_sig_aug(5,i);
    double nu_psy = x_sig_aug(6,i);

    Eigen::VectorXd x_vec = VectorXd(n_x_);
    // x_vec << px, py, v, psy, psy_dot;
    x_vec = x_sig_aug.col(i).head(5);

    //construct tranformation function
    double x_func, y_func;
    if(fabs(psy_dot) > 0.001)
    {
      x_func = (v/psy_dot)*( sin(psy+ psy_dot*delta_t) - sin(psy));
      y_func = (v/psy_dot)*(-cos(psy+ psy_dot*delta_t) + cos(psy));
    }
    else
    {
      x_func = v * cos(psy) * delta_t;
      y_func = v * sin(psy) * delta_t;
    }

    Eigen::VectorXd trans_func = VectorXd(n_x_);
    trans_func << x_func, y_func, 0, psy_dot*delta_t, 0.0;

    //noise vector
    Eigen::VectorXd noise = VectorXd(n_x_);
    noise << 0.5*nu_a*cos(psy)*delta_t*delta_t,
             0.5*nu_a*sin(psy)*delta_t*delta_t,
             delta_t*nu_a,
             0.5*delta_t*delta_t*nu_psy,
             delta_t*nu_psy;

    Xsig_pred_.col(i) = x_vec + trans_func + noise;
  }

  //predict the new mean value from sigma points
  Eigen::VectorXd pred_mean = VectorXd(n_x_);
  pred_mean.fill(0.0);
  x_.fill(0.0);
  for (int i = 0; i <=2*n_aug_; i++)
  {
    x_ += (weights_(i)*Xsig_pred_.col(i));
  }

  //predict the new covariance matrix from the sigma points
  Eigen::MatrixXd pred_covar = MatrixXd(n_x_, n_x_);
  pred_covar.fill(0.0);
  P_.fill(0.0);

  for (int i = 0; i <=2*n_aug_; i++)
  {
    Eigen::VectorXd mean_diff = Xsig_pred_.col(i) - x_;
    // std::cout << "diferenc " << mean_diff(3) << std::endl; 
    while (mean_diff(3)> M_PI) 
      mean_diff(3)-=2.*M_PI;
    while (mean_diff(3)<-M_PI) 
      mean_diff(3)+=2.*M_PI;

    P_ += (weights_(i)*mean_diff*mean_diff.transpose());
  }


  //now copy the predicted values as our new mean and covariance
  // x_ = pred_mean;
  // P_ = pred_covar;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
 
  //measurement size
  int n_z = 2;
  //create the new measurement sigma points based on the measurement function.
  Eigen::MatrixXd Zsigma = MatrixXd (n_z, 2*n_aug_+1);

  for(int i=0; i<=2*n_aug_; i++)
  {
    double x_pose = Xsig_pred_(0,i);
    double y_pose = Xsig_pred_(1,i);

    Zsigma(0,i) = x_pose;
    Zsigma(1,i) = y_pose;
  }


  //predicted the measurement state and covariance matrix
  Eigen::VectorXd z_mean_pred = VectorXd(n_z);
  z_mean_pred.fill(0.0);

  for (int i = 0; i<= 2*n_aug_; i++)
  {
    z_mean_pred += (weights_(i)*Zsigma.col(i));
  }

  Eigen::MatrixXd z_covar_pred = MatrixXd(n_z, n_z);
  z_covar_pred.fill(0.0);

  for (int i = 0; i<= 2*n_aug_; i++)
  {
    Eigen::VectorXd z_diff = Zsigma.col(i) - z_mean_pred;

    z_covar_pred += (weights_(i) * z_diff * z_diff.transpose());
  }

// std::cout << "diff: \n" << z_covar_pred << "\n =======+> " << std::endl;
  z_covar_pred += noise_lidar_;

  //calculate the cross-correlation matrix.
  Eigen::MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for  (int i = 0; i<= 2*n_aug_; i++)
  {

    //residual in measurement
    Eigen::VectorXd z_diff = Zsigma.col(i) - z_mean_pred;

    //difference in state prediction
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalizing the angle
    while (x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

    Tc += (weights_(i) * x_diff * z_diff.transpose());
  }

  //calculate kalman gain
  Eigen::MatrixXd gain = Tc * z_covar_pred.inverse();

  //now update the actual state vector and covariance matrix.
  Eigen::VectorXd input_meansurement = VectorXd(n_z);
  input_meansurement << meas_package.raw_measurements_(0),
                        meas_package.raw_measurements_(1);

  x_ += (gain * (input_meansurement - z_mean_pred));
  P_ -= (gain * z_covar_pred * gain.transpose());

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
 
    //measurement size
  int n_z = 3;
  //create the new measurement sigma points based on the measurement function.
  Eigen::MatrixXd Zsigma = MatrixXd (n_z, 2*n_aug_+1);

  for(int i=0; i<=2*n_aug_; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double theta = Xsig_pred_(3,i);

    double v1 = v*cos(theta);
    double v2 = v*sin(theta);

    Zsigma(0,i) = sqrt(px*px + py*py);
    Zsigma(1,i) = atan2(py, px);
    Zsigma(2,i) = (px*v1 + py*v2)/sqrt(px*px + py*py);
  }

  //predicted the measurement state and covariance matrix
  Eigen::VectorXd z_mean_pred = VectorXd(n_z);
  z_mean_pred.fill(0.0);

  for (int i = 0; i<= 2*n_aug_; i++)
  {
    z_mean_pred += (weights_(i)*Zsigma.col(i));
  }

  Eigen::MatrixXd z_covar_pred = MatrixXd(n_z, n_z);
  z_covar_pred.fill(0.0);

  for (int i = 0; i<= 2*n_aug_; i++)
  {
    Eigen::VectorXd z_diff = Zsigma.col(i) - z_mean_pred;

    // angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    z_covar_pred += (weights_(i)*z_diff*z_diff.transpose());
  }

  z_covar_pred += noise_radar_;

  //calculate the cross-correlation matrix.
  Eigen::MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for  (int i = 0; i<= 2*n_aug_; i++)
  {

    //residual in measurement
    Eigen::VectorXd z_diff = Zsigma.col(i) - z_mean_pred;

    // angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    //difference in state prediction
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalizing the angle
    while (x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

    Tc += (weights_(i) * x_diff * z_diff.transpose());
  }

  //calculate kalman gain
  Eigen::MatrixXd gain = Tc * z_covar_pred.inverse();

  //now update the actual state vector and covariance matrix.
  Eigen::VectorXd input_meansurement = VectorXd(n_z);
  input_meansurement << meas_package.raw_measurements_(0),
                        meas_package.raw_measurements_(1),
                        meas_package.raw_measurements_(2);

  Eigen::VectorXd residue = (input_meansurement - z_mean_pred);
  // angle normalization
  while (residue(1) > M_PI)  residue(1) -= 2. * M_PI;
  while (residue(1) < -M_PI) residue(1) += 2. * M_PI;

  x_ += gain * residue;
  P_ -= (gain * z_covar_pred * gain.transpose());

}