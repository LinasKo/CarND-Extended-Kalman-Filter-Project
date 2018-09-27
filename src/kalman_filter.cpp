#include "kalman_filter.h"

#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict()
{
    // DONE
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
    // DONE
    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    x_ += (K * y);
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
    // DONE: update the state by using Extended Kalman Filter equations
    MatrixXd Hj_ = tools_.CalculateJacobian(x_);
    if (Hj_.isZero(0))
    {
        return;
    }

    VectorXd y = z - tools_.CartesianToPolar(x_);
    MatrixXd S = Hj_ * P_ * Hj_.transpose() + R_;
    MatrixXd K = P_ * Hj_.transpose() * S.inverse();

    x_ += (K * y);
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I - K * Hj_) * P_;
}
