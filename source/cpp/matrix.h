#pragma once
#include<iostream>
using namespace std;
int sgn(double d);//符号函数long long factorial(long long n);//阶乘函数
class Matrix {
	//private:
public:
		int rows_num, cols_num;
		double** p;
		void initialize();//初始化矩阵

	//public:
		Matrix(int, int);
		Matrix();
		Matrix(int, int, double);//预配分空间
		Matrix(int m, int n, int k);//m*n,k阶单位阵,n=1
		Matrix(const Matrix& A);//复制构造函数
		~Matrix();//析构函数应当是虚函数，除非此类不用做基类
		Matrix& operator=(const Matrix&);//矩阵的复制
		Matrix& operator=(double*);//将数组的值传给矩阵
		Matrix& operator=(int*);//将数组的值传给矩阵
		Matrix& operator+=(const Matrix&);//矩阵的+=操作
		Matrix& operator-=(const Matrix&);//-=
		Matrix operator-(Matrix);
		Matrix& operator*=(const Matrix&);//*=
		Matrix operator*(const Matrix& m)const;
		Matrix operator*(double);
		Matrix operator/(double);
		Matrix operator+(const Matrix&);//+
		Matrix operator^(int);//重载^

		Matrix number_times(double m);//数乘
		static Matrix Solve(Matrix&, Matrix&);//求解线性方程组Ax=b
		static Matrix Solve_only(Matrix&, Matrix&);
		void Show() const;//矩阵显示
		void swapRows(int, int);
		double det();//求矩阵的行列式
		double Point(int i, int j) const;
		static Matrix inv(Matrix);//求矩阵的逆矩阵
		//static Matrix inv(const Matrix&);//求矩阵的逆矩阵
		static Matrix eye(int);//制造一个单位矩阵
		int row() const;
		int col() const;
		static Matrix T(const Matrix& m);//矩阵转置的实现,且不改变矩阵
		Matrix gaussianEliminate();//高斯消元法
		void gaussianEliminate(Matrix& b);//高斯消元法
		friend std::istream& operator>>(std::istream&, Matrix&);//重载>> 声明为友元函数访问行和列的值，实现矩阵的输入
		void householder(Matrix& b);//householder with b
		void householder();//
		Matrix* householder_m();
		double col_norm_2(int current_col, int k);//某一列第k个元素后的2范数
		void setpoint(int, int, double);//修改元素的值
		double norm2(int which_col);
		Matrix partial_Matrix(int row_begin, int row_end, int col_begin, int col_end);//得到矩阵的子矩阵
		static Matrix normalization(Matrix& A, Matrix& b);//正规化求解方程
		static Matrix augmented_solve(const Matrix& A, const Matrix& b);//增广矩阵求解
		static void givens_rotation(Matrix& A, Matrix& b);//givens rotation 
		static int find_nonzero(const Matrix& A, int which_col);//查找非零元素
		static int find_max(const Matrix& A, int which_col);
		static Matrix* gram_schmidt(Matrix& A, Matrix& b);//gram正交变换
		static Matrix* gram_schmidt(Matrix& A);//gram正交变换
		static Matrix* gram_improve(Matrix& A, Matrix& b);
		static Matrix* gram_improve(Matrix&);
		static Matrix* iteration_gram_schmidt(Matrix& A, Matrix& b);
		double infinite_norm(int);//某列的无穷范数
		Matrix exp_sum(int);//矩阵指数函数，int指定项数
		Matrix exp_gigenvector();//矩阵指数函数，int指定项数
		void gram_schmidt(Matrix& A, Matrix& Q, Matrix& R);//gram正交变换
		Matrix res(int,int);//余子式
		Matrix diag();
		Matrix combine_rows(Matrix& A, Matrix& B);
		Matrix combine_cols(Matrix& A, Matrix& B);

	};
	const double EPS = 0;

