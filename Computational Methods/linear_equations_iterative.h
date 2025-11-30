#pragma once
#include <vector>
#include <iostream>
#include <ranges> 
#include <algorithm>
using namespace std;
//线性方程组的迭代法
//Jacobi
//Gauss
//SOR
class iterative_method {
protected:
	vector<vector<double>> m_A;//迭代矩阵
	vector<double> m_B;
	vector<double> ans;
	double differ;//迭代次数
	unsigned n;
public:
	//构造函数
	iterative_method(const vector<vector<double>>& A, const vector<double>& b, double);
	//求解函数
	virtual void solve()=0;
	//打印输出结果
	void print() const;
	virtual ~iterative_method() {};
};
class Jacobi :public iterative_method {
public:
	//构造函数
	Jacobi(const vector<vector<double>>& A, const vector<double>& b, double);
	void solve() override;
};
class Gauss_Seidel :public iterative_method {
public:
	Gauss_Seidel(const vector<vector<double>>& A, const vector<double>& b, double);
	void solve() override;
};

//超松弛迭代
class SOR : public iterative_method {
private:
	double omega;  // 松弛因子
public:
	// 构造函数
	SOR(const vector<vector<double>>& A, const vector<double>& b, double diff, double w);
	void solve() override;
};
