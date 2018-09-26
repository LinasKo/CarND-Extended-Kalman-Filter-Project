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

    if (px == 0 and py == 0)
    {
        cout << "[Tools::CalculateJacobian] Error: px and py are both 0." << endl;
        return Hj;
    }

	//compute the Jacobian matrix
	float dist = sqrt(px * px + py * py);
	Hj << px / dist, py / dist, 0, 0,
	      -py / (dist * dist), px / (dist * dist), 0, 0,
	      py * (vx * py - vy * px) / (dist * dist * dist), px * (vy * px - vx * py) / (dist * dist * dist), px / dist, py / dist;

	return Hj;
}

VectorXd Tools::CartesianToPolar(const VectorXd& x_state)
{
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    VectorXd polar(3);
    polar <<
        sqrt(px * px + py * py),
        atan(py / px),
        (px * vx + py * vy) / sqrt(px * px + py * py);

    return polar;
}
