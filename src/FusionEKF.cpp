#include "FusionEKF.h"

#include <cmath>
#include <iostream>

#include "Eigen/Dense"
#include "tools.h"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
* Constructor.
*/
FusionEKF::FusionEKF()
{
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    /**
    * DONE:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
    */

    ekf_.x_ = VectorXd(4);

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;

    H_laser_  << 1, 0, 0, 0,
                 0, 1, 0, 0;
    ekf_.H_ = H_laser_;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
         	   0, 1, 0, 1,
         	   0, 0, 1, 0,
         	   0, 0, 0, 1;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_)
    {
        /**
        * DONE:
        * Initialize the state ekf_.x_ with the first measurement.
        * Create the covariance matrix.
        * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // Convert radar from polar to cartesian coordinates and initialize state
            float range = measurement_pack.raw_measurements_[0];
            float angle = measurement_pack.raw_measurements_[1];

            ekf_.x_ << range * cos(angle), range * sin(angle), 0, 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }

        // done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/

    /**
    * DONE:
    * Update the state transition matrix F according to the new elapsed time. - Time is measured in seconds.
    * Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	// dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    // Set up the F matrix
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // Set the process covariance
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4 / 4 * noise_ax_, 0, dt_3 / 2 * noise_ax_, 0,
               0, dt_4 / 4 * noise_ay_, 0, dt_3 / 2 * noise_ay_,
               dt_3 / 2 * noise_ax_, 0, dt_2 * noise_ax_, 0,
               0, dt_3 / 2 * noise_ay_, 0, dt_2 * noise_ay_;


    /*****************************************************************************
    *  Predict
    ****************************************************************************/

    ekf_.Predict();

    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
    * DONE:
    * Use the sensor type to perform the update step.
    * Update the state and covariance matrices.
    */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        cout << ">>>> RADAR" << endl;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else
    {
        cout << ">>>> LASER" << endl;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
