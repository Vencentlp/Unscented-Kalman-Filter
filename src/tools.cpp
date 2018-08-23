#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd RMSE(5);
	RMSE << 0, 0, 0, 0;
	if (estimations.size() != ground_truth.size())
	{
		cout << "Error - Calculate RMSE () - The groud truth and estimation vector must have same size";
		return RMSE;
	}
	for (int i = 0; i < estimations.size(); i++)
	{
		VectorXd Diff = estimations[i] - ground_truth[i];
		VectorXd residual = Diff.array() * Diff.array();
		RMSE = RMSE + residual;
	}
	RMSE = RMSE.array().sqrt();
	RMSE = RMSE / estimations.size();

	return RMSE;
}