#include "Unary_Equation.h"
#include <algorithm>
#include <cmath>
#include <iostream>
using namespace std;
Stationary_Points_Iterate::Stationary_Points_Iterate(double (*f)(double), double (*f2)(double)) :Function(f), g(f2) {}

double Stationary_Points_Iterate::solve(double x) const {
	double new_x = g(x);
	for(int i=0;i<1000;++i){
		cout << "µÚ" << i + 1 << "´Îµü´ú: " << new_x << endl;
		if (abs(new_x - x) < 1e-5)
			return new_x;
		x = new_x;
		new_x = g(x);
	}
	cout << "²»ÊÕÁ²£¡" << endl;
	return new_x;
}

double Newton::solve(double x)const {
	cout << "ÊäÈëÖØ¸ùÖØÊý£º" << std::endl;
	int m;
	std::cin >> m;
	return Nsolve(x, m);

}

double Newton::Nsolve(double x, int m)const {
	double new_x = x - 1.0 *m* func(x) / dfunc(x);
	for (int i = 0; i < 10; ++i) {
		cout << "µÚ" << i + 1 << "´Îµü´ú" << new_x << endl;
		if (abs(new_x - x) < 1e-5) {
			return new_x;
		}
		x = new_x;
		new_x = x - 1.0 *m* func(x) / dfunc(x);
	}
	cout << "²»ÊÕÁ²" << endl;
	return new_x;
}

double Steffensen::solve(double x)const {
	double y = g(x), z = g(y);
	double new_x = x - pow(y - x, 2) / (z - 2 * y + x);
	for (int i = 0; i < 10000; ++i) {
		if (abs(new_x - x) < 1e-5) {
			return new_x;
		}
		x = new_x;
		y = g(x), z = g(y);
		new_x = x - pow(y - x, 2) / (z - 2 * y + x);
		cout << "µÚ" << i + 1 << "´Îµü´ú: " << new_x << endl;
	}
	cout << "²»ÊÕÁ²!" << endl;
	return new_x;
}

