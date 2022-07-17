#pragma once
#include"complex.h"
#include"matrix.h"
class complex_matrix
{
public:
	int rows_num, cols_num;
	Complex** p;
	void initialize();//��ʼ������
public:
	/************���캯������������*/
	
	complex_matrix(int, int);
	complex_matrix();
	complex_matrix(const complex_matrix& );//���ƹ��캯��
	~complex_matrix();//��������


	/***********��������غͳ��ò���****************/

	Matrix realmatrix();//�õ�ʵ������
	Matrix imagmatrix();//�õ�ʵ������
	static complex_matrix eyes(int n);//��λ��
	complex_matrix operator=(const complex_matrix&);//����ĸ���
	complex_matrix& operator=(double*);//�������ֵ�������� 
	complex_matrix& operator=(Complex*);//�������ֵ�������� ������
	complex_matrix& operator+=(const complex_matrix&);//�����+=����
	complex_matrix& operator-=(const complex_matrix&);//-=
	complex_matrix operator-(const complex_matrix&);
	complex_matrix& operator*=(const complex_matrix&);//*=
	complex_matrix operator*(const complex_matrix& m);
	complex_matrix operator*(double);
	complex_matrix operator*(Complex);
	complex_matrix operator/(double);
	complex_matrix operator+(const complex_matrix&);//+
	complex_matrix operator^(int);//����^
	complex_matrix operator!()const;//����ת�õ�ʵ��,�Ҳ��ı����
	Complex operator()(int row, int col)const;//������ĳһ���ֵ
	complex_matrix col_vector(int which_col);//���������
	complex_matrix row_vector(int which_col);//���������
	complex_matrix inverse();//������
	complex_matrix pseudo_inverse()const;//α��
	complex_matrix& set_real_matrix(const Matrix&);//����ʵ������
	complex_matrix& set_imag_matrix(const Matrix&);//�����鲿����
	complex_matrix get_row(int j);
	complex_matrix get_col(int j);
	complex_matrix remove_row(int);
	complex_matrix remove_col(int);
	complex_matrix combine_rows(complex_matrix&,complex_matrix&);//�ϲ��������󣬴��з���
	complex_matrix combine_cols(complex_matrix&, complex_matrix&);//�ϲ��������󣬴��з���
	Complex cdot();//��ͬ�������ڻ�
	Complex cdot(complex_matrix);//��ͬ�������ڻ�
	static complex_matrix* gram_improve(complex_matrix& A);
	void exchange_col(int j1, int j2);//����
	void exchange_row(int i1, int i2);//����
	void Show() const;//������ʾ
	int row() const;
	int col() const;
	double norm2(int which_col);
	double norm2_vector();
};
	
