#pragma once
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;
class BaseCondition {
protected:
	double x0;
	double y0;
	double h;
	double xn;
	double (*dy)(double x, double y);
	vector<double> y_val;
public:
	BaseCondition(double, double, double (*d)(double x, double y),double h,double maxx);
	virtual void solve()= 0;
	vector<double> gety() const { return y_val; }
	virtual ~BaseCondition() {};
};

class Runge_Kutta4:public BaseCondition {
private:
	double k1(double x, double y)const { return dy(x, y); }
	double k2(double x, double y)const { return dy(x + 0.5 * h, y + 0.5 * h * k1(x, y)); }
	double k3(double x, double y)const { return dy(x + 0.5 * h, y + 0.5 * h * k2(x, y)); }
	double k4(double x, double y)const { return dy(x + h, y + h * k3(x, y)); }
public:
	Runge_Kutta4(double x, double y, double (*d)(double x, double y), double h, double maxx):BaseCondition(x,y,d,h,maxx){}
	void solve() override;
};