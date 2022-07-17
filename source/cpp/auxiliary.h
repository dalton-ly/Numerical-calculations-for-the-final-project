#pragma once
#include"complex_matrix.h"
#include<ctime>
#include<cstdlib>
/*******辅助函数头文件******/
int* generate_random_bits(int nT);
int* signal_to_bits(complex_matrix& c);
complex_matrix generate_H(int nR, int nT);
complex_matrix generate_signal(int nT);
complex_matrix generate_signal(int* ,int nT);
complex_matrix generate_noise(int nR);
//complex_matrix rand_matrix(const complex_matrix& A);
double gaussianrand(double=1.0 );
Complex quantization(Complex y);
complex_matrix quantization(complex_matrix);
double BitsErrorRate(complex_matrix, complex_matrix);
double partial(double (*f)(Matrix, Matrix), Matrix x, Matrix A, int i);//求f函数在第i个变量的偏导
Complex partial(Complex (*f)(Matrix, Matrix), Matrix x, Matrix A, int i);//求f函数在第i个变量的偏导
Matrix get_negetive_gradient(double (*f)(Matrix, Matrix,Matrix), Matrix x, Matrix A,Matrix);
complex_matrix get_negetive_gradient(Complex (*f)(complex_matrix, complex_matrix, complex_matrix,complex_matrix), complex_matrix x, complex_matrix A, complex_matrix,
	complex_matrix);//复数函数求梯度

complex_matrix gradient_complex(complex_matrix A, complex_matrix s, complex_matrix b);
complex_matrix Matrix_to_complex(Matrix);