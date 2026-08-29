#pragma once
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <utility>
using namespace std;
class Romberg {
private:
	double sup;
	double inf;
	function<double(double)> f;
	vector<vector<double>> R;
	double trapezoid(int k);
	int printn;
public:
	Romberg(double x,double y,function<double(double)> func):sup(x),inf(y),f(func),printn(0){
		if (x < y) {
			swap(sup, inf);
		}
	}
	double solve(int max_k = 10, double tolerance = 1e-4);
	void printTable() const;
};