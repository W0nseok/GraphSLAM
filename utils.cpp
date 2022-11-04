#include "utils.h"

Eigen::MatrixXd t2v(Eigen::MatrixXd t)
{
	Eigen::MatrixXd v(3, 1);
	v(0, 0) = t(0, 2);
	v(1, 0) = t(1, 2);
	v(2, 0) = atan2(t(1, 0), t(0, 0));

	return v;
}

Eigen::MatrixXd v2t(Eigen::MatrixXd v)
{
	Eigen::MatrixXd t(3, 3);
	double co = cos(v(2));
	double si = sin(v(2));
	t << co, -si, v(0, 0),
		si, co, v(1, 0),
		0, 0, 1;
	return t;
}