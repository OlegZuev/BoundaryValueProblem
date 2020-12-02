#pragma once
#include <map>
#include <ostream>

const int VARIANT = 5;
const double EPS_AUTO_STEP = 1E-5;
const double EPS_LIMIT = 1E-4;

double coef_A(double x);

double coef_B(double x);

double coef_C(double x);

double coef_F(double x);

double coef_Y_pr(double x, int derivative_number);

void runge_kutty_step(double z, double y, int i, double h, double& new_z, double& new_y);

void runge_kutty_method(double y0, double alpha0, double& h, int iter, double& max_delta,
                                            std::map<float, double>& z, std::map<float, double>& y, std::ostream& ostr);

void shooting_method(std::ostream& ostr);
