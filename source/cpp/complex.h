#pragma once
#include<iostream>
using namespace std;
class Complex
{
public:
	Complex()
	{ real =0; imag = 0; }
	Complex(double, double);

	//运算符重载
	Complex(const Complex&);
	Complex& operator=(const Complex&);
	Complex operator+=(Complex& a);
	friend Complex operator + (const Complex& c1,const Complex& c2);
	friend Complex operator + (double& d1, Complex& c2);
	friend Complex operator + (Complex& c1, double& d2);

	friend Complex operator - (const Complex& c1,const Complex& c2);
	friend Complex operator - (double& d1, Complex& c2);
	friend Complex operator - (Complex& c1, double& d2);

	friend Complex operator * (Complex& c1, Complex& c2);
	friend Complex operator * (double& d1, Complex& c2);
	friend Complex operator * (Complex& c1, double& d2);

	friend Complex operator / (const Complex& c1,const Complex& c2);
	friend Complex operator / (double& d1, Complex& c2);
	friend Complex operator / (Complex& c1, double& d2);

	//虚数单位用j表示
	friend ostream& operator << (ostream& out, Complex& c1);
	friend istream& operator >> (istream& in, Complex& c1);

	Complex operator!();

	//常用操作
	//求模值
	double modulus();	//为了单独访问实部

	double real, imag;
};

