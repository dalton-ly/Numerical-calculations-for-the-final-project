#pragma once
#include"complex.h"
#include"matrix.h"
class complex_matrix
{
public:
	int rows_num, cols_num;
	Complex** p;
	void initialize();//初始化矩阵
public:
	/************构造函数和析构函数*/
	
	complex_matrix(int, int);
	complex_matrix();
	complex_matrix(const complex_matrix& );//复制构造函数
	~complex_matrix();//析构函数


	/***********运算符重载和常用操作****************/

	Matrix realmatrix();//得到实部矩阵
	Matrix imagmatrix();//得到实部矩阵
	static complex_matrix eyes(int n);//单位阵
	complex_matrix operator=(const complex_matrix&);//矩阵的复制
	complex_matrix& operator=(double*);//将数组的值传给矩阵 
	complex_matrix& operator=(Complex*);//将数组的值传给矩阵 复数组
	complex_matrix& operator+=(const complex_matrix&);//矩阵的+=操作
	complex_matrix& operator-=(const complex_matrix&);//-=
	complex_matrix operator-(const complex_matrix&);
	complex_matrix& operator*=(const complex_matrix&);//*=
	complex_matrix operator*(const complex_matrix& m);
	complex_matrix operator*(double);
	complex_matrix operator*(Complex);
	complex_matrix operator/(double);
	complex_matrix operator+(const complex_matrix&);//+
	complex_matrix operator^(int);//重载^
	complex_matrix operator!()const;//矩阵转置的实现,且不改变矩阵
	Complex operator()(int row, int col)const;//返回在某一点的值
	complex_matrix col_vector(int which_col);//获得列向量
	complex_matrix row_vector(int which_col);//获得行向量
	complex_matrix inverse();//方阵逆
	complex_matrix pseudo_inverse()const;//伪逆
	complex_matrix& set_real_matrix(const Matrix&);//设置实部矩阵
	complex_matrix& set_imag_matrix(const Matrix&);//设置虚部矩阵
	complex_matrix get_row(int j);
	complex_matrix get_col(int j);
	complex_matrix remove_row(int);
	complex_matrix remove_col(int);
	complex_matrix combine_rows(complex_matrix&,complex_matrix&);//合并两个矩阵，从行方向
	complex_matrix combine_cols(complex_matrix&, complex_matrix&);//合并两个矩阵，从列方向
	Complex cdot();//相同向量的内积
	Complex cdot(complex_matrix);//不同向量的内积
	static complex_matrix* gram_improve(complex_matrix& A);
	void exchange_col(int j1, int j2);//换列
	void exchange_row(int i1, int i2);//换行
	void Show() const;//矩阵显示
	int row() const;
	int col() const;
	double norm2(int which_col);
	double norm2_vector();
};
	
