#pragma once
class Function {
protected:
	double (*func)(double x);
public:
	Function(double (*f)(double)) :func(f) {};
	virtual double solve(double x) const =0;
};
class Stationary_Points_Iterate :public Function {
protected:
	double (*g)(double x);
public:
	Stationary_Points_Iterate(double (*f)(double), double (*f2)(double));
	double solve(double x) const;
};
class Newton :public Function {
private:
	double (*dfunc)(double x);
	double Nsolve(double x, int m)const;
public:
	Newton(double (*f)(double), double (*df)(double)):Function(f),dfunc(df){}
	double solve(double x) const;
};
class Steffensen :public Stationary_Points_Iterate {
public:
	Steffensen(double (*f)(double), double (*f2)(double)) :Stationary_Points_Iterate(f, f2){}
	double solve(double x)const;
};
