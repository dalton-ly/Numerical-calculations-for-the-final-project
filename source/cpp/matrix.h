#pragma once
#include<iostream>
using namespace std;
int sgn(double d);//���ź���long long factorial(long long n);//�׳˺���
class Matrix {
	//private:
public:
		int rows_num, cols_num;
		double** p;
		void initialize();//��ʼ������

	//public:
		Matrix(int, int);
		Matrix();
		Matrix(int, int, double);//Ԥ��ֿռ�
		Matrix(int m, int n, int k);//m*n,k�׵�λ��,n=1
		Matrix(const Matrix& A);//���ƹ��캯��
		~Matrix();//��������Ӧ�����麯�������Ǵ��಻��������
		Matrix& operator=(const Matrix&);//����ĸ���
		Matrix& operator=(double*);//�������ֵ��������
		Matrix& operator=(int*);//�������ֵ��������
		Matrix& operator+=(const Matrix&);//�����+=����
		Matrix& operator-=(const Matrix&);//-=
		Matrix operator-(Matrix);
		Matrix& operator*=(const Matrix&);//*=
		Matrix operator*(const Matrix& m)const;
		Matrix operator*(double);
		Matrix operator/(double);
		Matrix operator+(const Matrix&);//+
		Matrix operator^(int);//����^

		Matrix number_times(double m);//����
		static Matrix Solve(Matrix&, Matrix&);//������Է�����Ax=b
		static Matrix Solve_only(Matrix&, Matrix&);
		void Show() const;//������ʾ
		void swapRows(int, int);
		double det();//����������ʽ
		double Point(int i, int j) const;
		static Matrix inv(Matrix);//�����������
		//static Matrix inv(const Matrix&);//�����������
		static Matrix eye(int);//����һ����λ����
		int row() const;
		int col() const;
		static Matrix T(const Matrix& m);//����ת�õ�ʵ��,�Ҳ��ı����
		Matrix gaussianEliminate();//��˹��Ԫ��
		void gaussianEliminate(Matrix& b);//��˹��Ԫ��
		friend std::istream& operator>>(std::istream&, Matrix&);//����>> ����Ϊ��Ԫ���������к��е�ֵ��ʵ�־��������
		void householder(Matrix& b);//householder with b
		void householder();//
		Matrix* householder_m();
		double col_norm_2(int current_col, int k);//ĳһ�е�k��Ԫ�غ��2����
		void setpoint(int, int, double);//�޸�Ԫ�ص�ֵ
		double norm2(int which_col);
		Matrix partial_Matrix(int row_begin, int row_end, int col_begin, int col_end);//�õ�������Ӿ���
		static Matrix normalization(Matrix& A, Matrix& b);//���滯��ⷽ��
		static Matrix augmented_solve(const Matrix& A, const Matrix& b);//����������
		static void givens_rotation(Matrix& A, Matrix& b);//givens rotation 
		static int find_nonzero(const Matrix& A, int which_col);//���ҷ���Ԫ��
		static int find_max(const Matrix& A, int which_col);
		static Matrix* gram_schmidt(Matrix& A, Matrix& b);//gram�����任
		static Matrix* gram_schmidt(Matrix& A);//gram�����任
		static Matrix* gram_improve(Matrix& A, Matrix& b);
		static Matrix* gram_improve(Matrix&);
		static Matrix* iteration_gram_schmidt(Matrix& A, Matrix& b);
		double infinite_norm(int);//ĳ�е������
		Matrix exp_sum(int);//����ָ��������intָ������
		Matrix exp_gigenvector();//����ָ��������intָ������
		void gram_schmidt(Matrix& A, Matrix& Q, Matrix& R);//gram�����任
		Matrix res(int,int);//����ʽ
		Matrix diag();
		Matrix combine_rows(Matrix& A, Matrix& B);
		Matrix combine_cols(Matrix& A, Matrix& B);

	};
	const double EPS = 0;

