#include "complex.h"
Complex::Complex(double re, double im)
{
	real = re;
	imag = im;
}
Complex::Complex(const Complex& a)
{
	real = a.real;
	imag = a.imag;
}
Complex& Complex::operator=(const Complex& a)
{
this->real = a.real;
this->imag = a.imag;
return *this;
}

Complex Complex::operator+=(Complex& a)
{
	this->real = this->real+a.real;
	this->imag = this->imag + a.imag;
	return *this;
}

Complex operator+(const Complex& c1, const Complex& c2)
{
	Complex c3;
	c3.real = c1.real + c2.real;
	c3.imag = c1.imag + c2.imag;
	return c3;
}
Complex operator + (double& d1, Complex& c2)
{
	Complex c3;
	c3.real = d1 + c2.real;
	c3.imag = c2.imag;
	return c3;
}
Complex operator + (Complex& c1, double& d2)
{
	Complex c3;
	c3.real = c1.real + d2;
	c3.imag = c1.imag;
	return c3;
}
//减号-
Complex operator - (const Complex& c1, const Complex& c2)
{
	Complex c3;
	c3.real = c1.real - c2.real;
	c3.imag = c1.imag - c2.imag;
	return c3;
}
Complex operator - (double& d1, Complex& c2)
{
	Complex c3;
	c3.real = d1 - c2.real;
	c3.imag = -c2.imag;
	return c3;
}
Complex operator - (Complex& c1, double& d2)
{
	Complex c3;
	c3.real = c1.real - d2;
	c3.imag = c1.imag;
	return c3;
}
//乘号*
Complex operator * (Complex& c1, Complex& c2)
{
	Complex c3;
	c3.real = c1.real * c2.real - c1.imag * c2.imag;
	c3.imag = c1.real * c2.imag + c1.imag * c2.real;
	return c3;
}
Complex operator * (double& d1, Complex& c2)
{
	Complex c3;
	c3.real = d1 * c2.real;
	c3.imag = d1 * c2.imag;
	return c3;
}
Complex operator * (Complex& c1, double& d2)
{
	Complex c3;
	c3.real = c1.real * d2;
	c3.imag = c1.imag * d2;
	return c3;
}
//除号/
Complex operator / (const Complex& c1,const Complex& c2)
{
	Complex c3;
	double temp = c2.real * c2.real + c2.imag * c2.imag;
	c3.real = (c1.real * c2.real + c1.imag * c2.imag) / temp;
	c3.imag = (c1.imag * c2.real - c1.real * c2.imag) / temp;
	return c3;
}
Complex operator / (double& d1, Complex& c2)
{
	Complex c3;
	double temp = c2.real * c2.real + c2.imag * c2.imag;
	c3.real = d1 * c2.real / temp;
	c3.imag = -d1 * c2.imag / temp;
	return c3;
}
Complex operator / (Complex& c1, double& d2)
{
	Complex c3(c1.real / d2, c1.imag / d2);
	return c3;
}

//输出<<,输入>>
ostream& operator << (ostream& out, Complex& c1)
{
	if (c1.imag >=0)
		out << c1.real << "+" << c1.imag << "j";
	else if (c1.imag < 0)
		out << c1.real << c1.imag << "j";
	return out;
}
istream& operator >> (istream& in, Complex& c1)	//先输入实部，后输入虚部
{
	in >> c1.real >> c1.imag;
	return in;
}
//共轭!
Complex Complex::operator !()
{
	Complex c1;
	c1.real = this->real;
	c1.imag = -(this->imag);
	return c1;
}

//求模值
double Complex::modulus()
{
	return sqrt(real * real + imag * imag);
}
