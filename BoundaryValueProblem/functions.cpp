#include "functions.h"
#include <map>
#include <ostream>
#include <iomanip>
#include <stdexcept>

double coef_A(double x) {
	switch (VARIANT) {
	case 0:
		return 0;
	case 5:
		return 40 * (x + 1);
	case 7:
		return 30 * (x - 0.5);
	}

	throw std::runtime_error("Incorrect variant");
}

double coef_B(double x) {
	switch (VARIANT) {
	case 0:
		return 0;
	case 5:
		return pow(x, 2) + 2;
	case 7:
		return pow(x, 2) + 1;
	}

	throw std::runtime_error("Incorrect variant");
}

double coef_C(double x) {
	switch (VARIANT) {
	case 0:
		return 0;
	case 5:
		return x + 1;
	case 7:
		return 2 * x - 1;
	}

	throw std::runtime_error("Incorrect variant");
}

double coef_F(double x) {
	return coef_Y_pr(x, 2) + coef_A(x) * coef_Y_pr(x, 1) - coef_B(x) * coef_Y_pr(x, 0) + coef_C(x) * sin(
		coef_Y_pr(x, 0));
}

double coef_Y_pr(double x, int derivative_number) {
	switch (derivative_number) {
	case 0:
		return 1.0 + x + 10 * log(VARIANT + 1) * pow(x, 3) * pow(1 - x, 3);
	case 1:
		return 1.0 - 30 * log(VARIANT + 1) * pow(x, 2) * pow(1 - x, 2) * (2 * x - 1);
	case 2:
		return -60.0 * x * log(VARIANT + 1) * (5 * pow(x, 3) - 10 * pow(x, 2) + 6 * x - 1);
	}

	throw std::runtime_error("Can`t compute this derivative of Ypr!");
}

void runge_kutty_step(double z, double y, int i, double h, double& new_z, double& new_y) {
	double k_y[4];
	double k_z[4];
	k_y[0] = h * z;
	k_z[0] = h * (-coef_A(i * h) * z + coef_B(i * h) * y - coef_C(i * h) * sin(y) + coef_F(i * h));
	k_y[1] = h * (z + k_z[0] / 2.0);
	k_z[1] = h * (-coef_A(i * h + h / 2.0) * (z + k_z[0] / 2.0) + coef_B(i * h + h / 2.0) * (y + k_y[0] / 2.0) -
		coef_C(i * h + h / 2.0) * sin(y + k_y[0] / 2.0) + coef_F(i * h + h / 2.0));
	k_y[2] = h * (z + k_z[1] / 2.0);
	k_z[2] = h * (-coef_A(i * h + h / 2.0) * (z + k_z[1] / 2.0) + coef_B(i * h + h / 2.0) * (y + k_y[1] / 2.0) -
		coef_C(i * h + h / 2.0) * sin(y + k_y[1] / 2.0) + coef_F(i * h + h / 2.0));
	k_y[3] = h * (z + k_z[2]);
	k_z[3] = h * (-coef_A(i * h + h) * (z + k_z[2]) + coef_B(i * h + h) * (y + k_y[2]) - coef_C(i * h + h) *
		sin(y + k_y[2]) + coef_F(i * h + h));

	new_y = y + 1.0 / 6.0 * (k_y[0] + 2.0 * k_y[1] + 2.0 * k_y[2] + k_y[3]);
	new_z = z + 1.0 / 6.0 * (k_z[0] + 2.0 * k_z[1] + 2.0 * k_z[2] + k_z[3]);
}

void runge_kutty_method(double y0, double alpha0, double& h, int iter, std::map<float, double>& z, std::map<float, double>& y, std::ostream& ostr) {
	z.clear();
	y.clear();
	double delta = -1.0;
	int i = 0;
	y[0] = y0;
	z[0] = alpha0;

	do {
		double z_j, y_j, z_half, y_half, z_asterisk, y_asterisk;

		runge_kutty_step(z[h * i], y[h * i], i, h, z_j, y_j);
		runge_kutty_step(z[h * i], y[h * i], i, h / 2.0, z_half, y_half);
		runge_kutty_step(z_half, y_half, i, h / 2.0, z_asterisk, y_asterisk);

		if (abs(y_j - y_asterisk) >= EPS_AUTO_STEP) {
			h /= 2.0;
			i = 0;
		} else {
			i++;
			y[i * h] = y_j;
			z[i * h] = z_j;
		}
	} while (i * h < 1.0);

	double p = 0.00625;
	for (double q = 0; q <= 1; q += p) {
		if (delta < abs(y[q] - coef_Y_pr(q, 0))) {
			delta = abs(y[q] - coef_Y_pr(q, 0));
		}

		if (abs(q - 0.96250) < 0.0001) {
			p /= 2;
		}
	}

	ostr << std::setw(4) << iter << "|" << std::setprecision(5) << std::fixed << std::setw(9) << z[0] << "|" <<
		std::setw(9) << y[1] << "|" << std::setw(19) << std::setprecision(11) << std::scientific << delta << "|" <<
		std::endl;
}


void shooting_method(std::ostream& ostr)
{
	double y0 = 1.0, h = 0.1, a0 = 0, a1 = 3.0;
	std::map<float, double> z, y;
	int i = 1;

	ostr << std::setw(4) << std::setprecision(0) << "Itr|" << std::setw(9) << "z(0)|"
		<< std::setw(9) << "y(0)|" << std::setw(20) << "Delta" << std::endl;

	runge_kutty_method(y0, a0, h, i++, z, y, ostr);
	runge_kutty_method(y0, a1, h, i++, z, y, ostr);

	do {
		double a2 = (a0 + a1) / 2;
		runge_kutty_method(y0, a2, h, i, z, y, ostr);
		if (y[1.0] < 2.0) {
			a0 = a2;
		} else {
			a1 = a2;
		}

		i++;
	} while (abs(y[1.0] - 2.0) > EPS_LIMIT);

	ostr << std::endl << std::setw(9) << "x|" << std::setw(9) << "y(x)|" << std::setw(9) << "Ypr|" << std::setw(9) << "z(x)|" << std::setw(19) << "Delta|" << std::endl;

	double max_delta = 0;
	double p = 0.00625;
	for (double q = 0; q <= 1; q += p)
	{
		double error = abs(y[q] - coef_Y_pr(q, 0));

		if (error > max_delta) {
			max_delta = error;
		}
		
		ostr << std::setw(9) << std::fixed << std::setprecision(5) << q << "|" << std::setw(9) << y[q]
			<< "|" << std::setw(9) << coef_Y_pr(q, 0) << "|" << std::setw(9) << z[q] << "|" << std::setw(19);

		if (max_delta < 0.0001) {
			ostr << std::setprecision(11) << std::scientific << error << "|" << std::endl;
		}
		else {
			ostr << std::setprecision(15) << std::fixed << error << "|" << std::endl;
		}
		
		if (abs(q - 0.96250) < 0.0001) {
			p /= 2;
		}
	}

	ostr << std::endl << "The absolute value of the maximum error: " << std::fixed << max_delta << std::endl;
}