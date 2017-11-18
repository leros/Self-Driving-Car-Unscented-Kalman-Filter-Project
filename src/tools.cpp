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
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  //check the validity of the inputs
  int esti_size = estimations.size();
  if((esti_size == 0) || (esti_size != ground_truth.size())){
	return rmse;
  }

  for(int i=0; i < esti_size; ++i){
	VectorXd diff = estimations[i] - ground_truth[i];
	diff = diff.array() * diff.array();
	rmse += diff;
  }

  rmse = rmse/esti_size;
  rmse = rmse.array().sqrt();

  return rmse;
}

