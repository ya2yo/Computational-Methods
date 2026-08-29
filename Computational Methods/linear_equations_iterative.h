#pragma once
#include <vector>
#include <iostream>
#include <ranges> 
#include <algorithm>
using namespace std;
//ÏßÐÔ·½³Ì×éµÄµü´ú·¨
//Jacobi
//Gauss
//SOR
class iterative_method {
protected:
	vector<vector<double>> m_A;//µü´ú¾ØÕó
	vector<double> m_B;
	vector<double> ans;
	double differ;//µü´ú´ÎÊý
	unsigned n;
public:
	//¹¹Ôìº¯Êý
	iterative_method(const vector<vector<double>>& A, const vector<double>& b, double);
	//Çó½âº¯Êý
	virtual void solve()=0;
	//´òÓ¡Êä³ö½á¹û
	void print() const;
	virtual ~iterative_method() {};
};
class Jacobi :public iterative_method {
public:
	//¹¹Ôìº¯Êý
	Jacobi(const vector<vector<double>>& A, const vector<double>& b, double);
	void solve() override;
};
class Gauss_Seidel :public iterative_method {
public:
	Gauss_Seidel(const vector<vector<double>>& A, const vector<double>& b, double);
	void solve() override;
};

//³¬ËÉ³Úµü´ú
class SOR : public iterative_method {
private:
	double omega;  // ËÉ³ÚÒò×Ó
public:
	// ¹¹Ôìº¯Êý
	SOR(const vector<vector<double>>& A, const vector<double>& b, double diff, double w);
	void solve() override;
};
