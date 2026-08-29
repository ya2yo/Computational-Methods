#include "Ordinary_Differential_Equations.h"

BaseCondition::BaseCondition(double x, double y, double(*f)(double x, double y), double h, double maxn)
    : x0(x), y0(y), h(h), xn(maxn), dy(f) {
	int n = round((xn - x0) / h);
	y_val.resize(n + 1);
	y_val[0] = y0;
}

void Runge_Kutta4::solve(){
	int n = y_val.size()-1;
	for (int i = 1; i <= n; ++i) {
		double x_prev = x0 + 1.0 * (i-1) * h,y_prev = y_val[i - 1];
		y_val[i] = y_prev + h / 6.0 * (k1(x_prev, y_prev) + 2 * k2(x_prev, y_prev) + 2 * k3(x_prev, y_prev) + k4(x_prev, y_prev));
	}
	for (int i = 0; i <= n; ++i) {
		double current_x = x0 + 1.0 * i * h;
		cout << "x[" << i << "] = " << current_x << ", y[" << i << "] = " << y_val[i] << endl;
	}
}
