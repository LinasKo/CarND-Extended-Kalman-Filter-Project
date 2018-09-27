#include "tools.h"

#include <cmath>
#include <iostream>


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}
Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
    VectorXd rmse(VectorXd::Zero(4));

    if (estimations.size() == 0)
    {
        std::cout << "[Tools::CalculateRMSE] Error: estimations length is 0" << std::endl;
        return rmse;
    }
    if (estimations.size() != ground_truth.size())
    {
        std::cout << "[Tools::CalculateRMSE] Error: estimations length differs from ground truth data length." << std::endl;
        return rmse;
    }

	for (size_t i = 0; i < estimations.size(); ++i)
    {
        VectorXd err = estimations[i] - ground_truth[i];
        VectorXd errSquared = err.array() * err.array();
		rmse += errSquared;
	}
    rmse /= estimations.size();
	rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
    MatrixXd Hj(3, 4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

    float sq_dist = px * px + py * py;
    float dist = sqrt(sq_dist);
    if (sq_dist < 0.0001)
    {
        cout << "[Tools::CalculateJacobian] Error: (px^2 + py^2) is too close to zero." << endl;
        Hj.setZero();
        return Hj;
    }

	//compute the Jacobian matrix
	Hj << px / dist, py / dist, 0, 0,
	      -py / sq_dist, px / sq_dist, 0, 0,
	      py * (vx * py - vy * px) / (dist * sq_dist), px * (vy * px - vx * py) / (sq_dist * dist), px / dist, py / dist;

	return Hj;
}

VectorXd Tools::CartesianToPolar(const VectorXd& x_state)
{
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    VectorXd polar(3);
    float dist = sqrt(px * px + py * py);
    float angle = atan2(py, px);
    float angular_vel = (px * vx + py * vy) / dist;
    polar << dist, angle, angular_vel;

    return polar;
}
