#include "matrix.h"
#include <cmath>
#include <stdlib.h>
#include<iomanip>
#include<algorithm>
void Matrix::initialize() //��ʼ�������С
{
	p = new double* [rows_num];//����rows_num��ָ��
	for (int i = 0; i < rows_num; ++i)
	{
		p[i] = new double[cols_num];//Ϊp[i]���ж�̬�ڴ���䣬��СΪcols
	}
}
Matrix::Matrix()
{
	rows_num = 1;
	cols_num = 1;
	initialize();

}
//����һ��ȫ0����
Matrix::Matrix(int rows, int cols)
{
	rows_num = rows;
	cols_num = cols;
	initialize();
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			p[i][j] = 0;
		}
	}
}
//����һ��ֵȫ��Ϊvalue�ľ���
Matrix::Matrix(int rows, int cols, double value)
{
	rows_num = rows;
	cols_num = cols;
	initialize();
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			p[i][j] = value;
		}
	}
}

//��������
Matrix::~Matrix() {
	for (int i = 0; i < rows_num; ++i)
	{
		delete[] p[i];
	}
	delete[] p;
}

//ʵ�־���ĸ���
Matrix& Matrix::operator=(const Matrix& m)
{
	if (this == &m)
	{
		return *this;
	}
	if (rows_num != m.rows_num || cols_num != m.cols_num)
	{
		//cout << "ά�Ȳ���ͬ" << endl;
		for (int i = 0; i < rows_num; ++i) {
			delete[] p[i];
		}
		delete[] p;

		rows_num = m.rows_num;
		cols_num = m.cols_num;
		initialize();
	}

	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			p[i][j] = m.p[i][j];
		}
	}
	return *this;
}
//�������ֵ���ݸ�����
Matrix& Matrix::operator=(double* a) 
{
	for (int i = 0; i < rows_num; i++) 
	{
		for (int j = 0; j < cols_num; j++) 
		{
			p[i][j] = *(a + i * cols_num + j);
		}
	}
	return *this;
}
Matrix& Matrix::operator=(int* a)
{
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			p[i][j] = double(*(a + i * cols_num + j));
		}
	}
	return *this;
}
//+=����
Matrix& Matrix::operator+=(const Matrix& m)
{
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] += m.p[i][j];
		}
	}
	return *this;
}
//ʵ��-=
Matrix& Matrix::operator-=(const Matrix& m)
{
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] -= m.p[i][j];
		}
	}
	return *this;
}
//ʵ��*=
Matrix& Matrix::operator*=(const Matrix& m)
{
	Matrix temp(rows_num, m.cols_num);//��C=AB,�����C���������ھ���A��������C����������B��������
	for (int i = 0; i < temp.rows_num; i++) {
		for (int j = 0; j < temp.cols_num; j++) {
			for (int k = 0; k < cols_num; k++) {
				temp.p[i][j] += (p[i][k] * m.p[k][j]);
			}
		}
	}
	*this = temp;
	return *this;
}
//ʵ�־���ĳ˷�
Matrix Matrix::operator*(const Matrix& m)const
{
	if (cols_num != m.rows_num)
	{
		cout << "ά�Ȳ���" << endl;
	}
	Matrix ba_M(rows_num, m.cols_num, 0.0);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < m.cols_num; j++)
		{
			for (int k = 0; k < cols_num; k++)
			{
				ba_M.p[i][j] += (p[i][k] * m.p[k][j]);
			}
		}
	}
	return ba_M;
}

Matrix Matrix::Solve(Matrix& A, Matrix& b)
{

	//��˹��ȥ��ʵ��Ax=b�ķ������
	for (int i = 0; i < A.rows_num; i++)
	{

		if (abs(A.p[i][i]) < EPS) //��Ҫ����ѡ��Ԫ
		{
			cout << "��ԪΪ0��ѡ��Ԫ" << endl;

		}
		for (int j = i + 1; j < A.rows_num; j++)
		{
			for (int k = i + 1; k < A.cols_num; k++)
			{
				A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
				if (abs(A.p[j][k]) < EPS)
					A.p[j][k] = 0;
			}
			b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
			if (abs(A.p[j][0]) < EPS)
				A.p[j][0] = 0;
			A.p[j][i] = 0;
		}

	}
	//A.Show();

	// �������
	Matrix x(b.rows_num, 1);
	x.p[x.rows_num - 1][0] = b.p[x.rows_num - 1][0] / A.p[x.rows_num - 1][x.rows_num - 1];
	if (abs(x.p[x.rows_num - 1][0]) < EPS)
		x.p[x.rows_num - 1][0] = 0;
	for (int i = x.rows_num - 2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < x.rows_num; j++)
		{
			sum += A.p[i][j] * x.p[j][0];
		}
		x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
		if (abs(x.p[i][0]) < EPS)
			x.p[i][0] = 0;
	}

	return x;
}


/*Matrix Matrix::Solve( Matrix& A, Matrix& b)
{

	//��˹��ȥ��ʵ��Ax=b�ķ������
	for (int i = 0; i < A.rows_num; i++)
	{
		A.Show();
		if (abs(A.p[i][i])<EPS) //��Ҫ����ѡ��Ԫ
		{

			cout << "����������" << endl;
		}
		for (int j = i + 1; j < A.rows_num; j++)
		{
			for (int k = i + 1; k < A.cols_num; k++)
			{
				A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
				if (abs(A.p[j][k]) < EPS)
					A.p[j][k] = 0;
			}
			b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
			if (abs(A.p[j][0]) < EPS)
				A.p[j][0] = 0;
			A.p[j][i] = 0;
		}
	}
	//A.Show();

	// �������
	Matrix x(b.rows_num, 1);
	x.p[x.rows_num - 1][0] = b.p[x.rows_num - 1][0] / A.p[x.rows_num - 1][x.rows_num - 1];
	if (abs(x.p[x.rows_num - 1][0]) <EPS)
		x.p[x.rows_num - 1][0] = 0;
	for (int i = x.rows_num - 2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < x.rows_num; j++)
		{
			sum += A.p[i][j] * x.p[j][0];
		}
		x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
		if (abs(x.p[i][0]) <EPS)
			x.p[i][0] = 0;
	}

	return x;
}
*/

//������ʾ
void Matrix::Show() const {
	//cout << rows_num <<" "<<cols_num<< endl;//��ʾ���������������
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			cout << setw(7) << setprecision(10) << fixed << " " << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
//ʵ���б任
void Matrix::swapRows(int a, int b)
{
	a--;
	b--;
	double* temp = p[a];
	p[a] = p[b];
	p[b] = temp;
}
//�����������ʽ��ֵ
/*double Matrix::det() {
	//Ϊ��������ʽ��һ������
	double** back_up;
	back_up = new double* [rows_num];
	for (int i = 0; i < rows_num; i++) {
		back_up[i] = new double[cols_num];
	}
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			back_up[i][j] = p[i][j];
		}
	}
	if (rows_num != cols_num) {
		std::abort();//ֻ�з�����ܼ�������ʽ����������ж�ǿ��ֹͣ����
	}
	double ans = 1;
	for (int i = 0; i < rows_num; i++)
	{
		//ͨ���б仯����ʽ��ʹ�þ���Խ����ϵ���Ԫ�ز�Ϊ0
		if (abs(p[i][i]) <= EPS) {
			bool flag = false;
			for (int j = 0; (j < cols_num) && (!flag); j++) {
				//�������һ���Խ����ϵ�Ԫ�ؽӽ���0���ܹ�ͨ���б任ʹ�þ���Խ����ϵ�Ԫ�ز�Ϊ0
				if ((abs(p[i][j]) > EPS) && (abs(p[j][i]) > EPS)) {
					flag = true;
					//�Ծ�������б任
					double temp;
					for (int k = 0; k < cols_num; k++) {
						temp = p[i][k];
						p[i][k] = p[j][k];
						p[j][k] = temp;
					}
				}
			}
			if (flag)
				return 0;
		}
	}
	for (int i = 0; i < rows_num; i++) {
		for (int j = i + 1; j < rows_num; j++) {
			for (int k = i + 1; k < cols_num; k++) {
				p[j][k] -= p[i][k] * (p[j][i] * p[i][i]);
			}
		}
	}
	for (int i = 0; i < rows_num; i++) {
		ans *= p[i][i];
	}
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] = back_up[i][j];
		}
	}
	return ans;
}
*/
//���ؾ����i�е�j�е���
double Matrix::Point(int i, int j) const
{
	return this->p[i][j];
}
//�����������
///*
Matrix Matrix::inv(Matrix A)
{
	if (A.rows_num != A.cols_num)
	{
		std::cout << "ֻ�з������������" << std::endl;
		std::abort();//ֻ�з������������
	}
	
	double temp;
	Matrix A_B = Matrix(A.rows_num, A.cols_num);
	A_B = A;//Ϊ����A��һ������
	Matrix B = eye(A.rows_num);
	//��С��EPS����ȫ����0
	for (int i = 0; i < A.rows_num; i++)
	{
		for (int j = 0; j < A.cols_num; j++)
		{
			if (abs(A.p[i][j]) <= EPS)
			{
				A.p[i][j] = 0;
			}
		}
	}
	//ѡ����Ҫ����������ѡ��Ԫ
	for (int i = 0; i < A.rows_num; i++)
	{
		if (abs(A.p[i][i]) <= EPS)
		{
			//bool flag = false;
			for (int j = 0; (j < A.rows_num) /*&& (!flag)*/; j++)
			{
				if ((abs(A.p[i][j]) > EPS) && (abs(A.p[j][i]) > EPS))
				{
					//flag = true;
					for (int k = 0; k < A.cols_num; k++)
					{
						temp = A.p[i][k];
						A.p[i][k] = A.p[j][k];
						A.p[j][k] = temp;
						temp = B.p[i][k];
						B.p[i][k] = B.p[j][k];
						B.p[j][k] = temp;
					}
				}
			}
			/*if (!flag)
			{
					std::cout << "����󲻴���\n";
					std::abort();
			}*/
		}
	}
	//ͨ�������б任��A��Ϊ�����Ǿ���
	double temp_rate;
	for (int i = 0; i < A.rows_num; i++) {
		for (int j = i + 1; j < A.rows_num; j++) {
			temp_rate = A.p[j][i] / A.p[i][i];
			for (int k = 0; k < A.cols_num; k++) {
				A.p[j][k] -= A.p[i][k] * temp_rate;
				B.p[j][k] -= B.p[i][k] * temp_rate;
			}
			A.p[j][i] = 0;
		}
	}
	//ʹ�Խ�Ԫ�ؾ�Ϊ1
	for (int i = 0; i < A.rows_num; i++) {
		temp = A.p[i][i];
		for (int j = 0; j < A.cols_num; j++) {
			A.p[i][j] /= temp;
			B.p[i][j] /= temp;
		}
	}
	//���Ѿ���Ϊ�����Ǿ����A����Ϊ��λ����
	for (int i = A.rows_num - 1; i >= 1; i--) {
		for (int j = i - 1; j >= 0; j--) {
			temp = A.p[j][i];
			for (int k = 0; k < A.cols_num; k++) {
				A.p[j][k] -= A.p[i][k] * temp;
				B.p[j][k] -= B.p[i][k] * temp;
			}
		}
	}
	// std::cout << "�㷨�ɿ��Լ�⣬���ɿ��������λ����" << std::endl;
	/*for (int i = 0; i < A.rows_num; i++) {
		for (int j = 0; j < A.cols_num; j++) {
			printf("%7.4lf\t\t", A.p[i][j]);
		}
		cout << endl;
	}*/
	A = A_B;//��ԭA
	return B;//���ظþ���������
}
//Matrix Matrix::inv( Matrix A)
//{
//	if (A.det() ==0)
//	{
//		cout << "����ʽΪ0������" << endl;
//	}
//	else
//	{
//		Matrix invmatrix(A.rows_num, A.cols_num);
//		for (int i = 0; i < A.rows_num; i++)
//		{
//			for (int j = 0; j < A.cols_num; j++)
//			{
//				if((i+j)%2==0)
//				invmatrix.setpoint(i, j, A.res(j, i).det());
//				else
//				invmatrix.setpoint(i, j, -A.res(j, i).det());
//			}
//		}
//		return invmatrix/A.det();
//	}
//}
 double Matrix::det()
{
	if(rows_num!=cols_num)
    {
		cout << "�Ƿ���������" << endl;
		abort();
    }
	else
	{
		if (rows_num == 1)
		{
			return this->Point(0, 0);
		}
		else
		{
			double det = 0;
			for (int i = 0; i < cols_num; i++)//iΪ��
			{
				Matrix temp(rows_num-1, cols_num-1);
				double a = this->Point(0, i);
				int k = 1;
				for (int temp_row = 0; temp_row < temp.rows_num; temp_row++)
				{
					
					for (int temp_col = 0; temp_col < temp.cols_num; temp_col++)
					{
						if (temp_col < i)
						{
							temp.setpoint(temp_row, temp_col, this->Point(k, temp_col));
						}
						else
						{
							temp.setpoint(temp_row, temp_col, this->Point(k, temp_col + 1));
						}
						
					}		
					k++;
				}
				
				if( i%2==0)
					det +=a*temp.det();
				else
					det-= a*temp.det();
			}
			return det;
		}
	}
}
//����һ����λ��
Matrix Matrix::eye(int n) {
	Matrix A(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				A.p[i][j] = 1.0;
			}
			else {
				A.p[i][j] = 0;
			}
		}
	}
	return A;
}
//����һ��m�е�����������k��Ԫ��Ϊ1
Matrix::Matrix(int m, int n, int k)
{
	rows_num = m;
	cols_num = 1;
	initialize();
	for (int i = 0; i < rows_num; i++)
	{
		p[i][0] = 0;

	}
	p[k][0] = 1;
}
//��ȡ����������
int Matrix::row() const {
	return rows_num;
}
int Matrix::col() const {
	return cols_num;
}
//ʵ�־����ת��
Matrix Matrix::T(const Matrix& m)
{
	int col_size = m.col();
	int row_size = m.row();
	Matrix mt(col_size, row_size);
	for (int i = 0; i < row_size; i++) {
		for (int j = 0; j < col_size; j++) {
			mt.p[j][i] = m.p[i][j];
		}
	}
	return mt;
}
//��˹��Ԫ��
Matrix Matrix::gaussianEliminate()
{
	Matrix Ab(rows_num, cols_num);
	Ab = *this;
	Ab.Show();
	Matrix m(rows_num, rows_num);

	for (int k = 0; k < rows_num - 1; k++)
	{
		m = Matrix::eye(rows_num);
		int maxrow = Matrix::find_max(Ab, k);
		if (maxrow != k)
		{
			Ab.swapRows(k, maxrow);
		}
		if (Ab.Point(k, k) == 0)
		{
			continue;
		}
		for (int i = k + 1; i < rows_num; i++)
		{
			m.setpoint(i, k, Ab.Point(i, k) / Ab.Point(k, k));
		}
		Ab = m * Ab;

	}

	return Ab;




	/*
	Matrix Ab(*this);
	int rows = Ab.rows_num;
	int cols = Ab.cols_num;
	int Acols = cols - 1;

	int i = 0; //������
	int j = 0; //������
	while (i < rows)
	{
		bool flag = false;
		while (j < Acols && !flag)
		{
			if (Ab.p[i][j] != 0)
			{
				flag = true;
			}
			else
			{
				int max_row = i;
				double max_val = 0;
				for (int k = i + 1; k < rows; ++k)
				{
					double cur_abs = Ab.p[k][j] >= 0 ? Ab.p[k][j] : -1 * Ab.p[k][j];
					if (cur_abs > max_val)
					{
						max_row = k;
						max_val = cur_abs;
					}
				}
				if (max_row != i)
				{
					Ab.swapRows(max_row, i);
					flag = true;
				}
				else
				{
					j++;
				}
			}
		}
		if (flag)
		{
			for (int t = i + 1; t < rows; t++)
			{
				for (int s = j + 1; s < cols; s++)
				{
					Ab.p[t][s] = Ab.p[t][s] - Ab.p[i][s] * (Ab.p[t][j] / Ab.p[i][j]);
					if (abs(Ab.p[t][s]) <= EPS)
						Ab.p[t][s] = 0;
				}
				Ab.p[t][j] = 0;
			}
		}
		i++;
		j++;
	}
	return Ab;
	*/
}
//ʵ�־��������
istream& operator>>(istream& is, Matrix& m)
{
	for (int i = 0; i < m.rows_num; i++) {
		for (int j = 0; j < m.cols_num; j++) {
			is >> m.p[i][j];
		}
	}
	return is;
}
//householder��b����
void Matrix::householder(Matrix& b)
{
	Matrix Vk(rows_num, 1);
	Matrix aj(rows_num, 1);
	for (int k = 0; k < cols_num; k++)
	{
		double alphak = -sgn(this->Point(k, k)) * (this->col_norm_2(k, k));
		//��Vk��ֵ
		for (int i = 0; i < rows_num; i++)
		{
			if (i < k)
			{
				Vk.setpoint(i, 0, 0);
			}
			else if (i == k)
			{
				Vk.setpoint(i, 0, (this->Point(k, k)) - alphak);
			}
			else if (i > k)
			{
				Vk.setpoint(i, 0, this->Point(i, k));
			}
		}
		double betaj = Vk.norm2(0);
		if (betaj == 0)
		{
			continue;
		}
		bool flag = false;//�ж��Ƿ��b��Ԫ
		for (int j = k; j < cols_num; j++)//������Ԫ
		{
			//�Ӿ����еõ�������
			for (int ii = 0; ii < rows_num; ii++)
			{
				aj.setpoint(ii, 0, this->Point(ii, j));
			}
			Matrix gamaj(1, 1, 0);
			Matrix gamaj_tob(1, 1, 0);
			Matrix Vkt(Vk.col(), Vk.row());
			Vkt = Matrix::T(Vk);
			gamaj = Vkt * aj;
			gamaj_tob = Vkt * b;
			//aj��Ԫ
			for (int jj = 0; jj < rows_num; jj++)
			{
				this->setpoint(jj, j, this->Point(jj, j) - 2 * (gamaj.Point(0, 0) / betaj) * Vk.Point(jj, 0));

			}
			if (flag == false)//��b��Ԫ
			{

				for (int jj = 0; jj < rows_num; jj++)
				{
					b.setpoint(jj, 0, b.Point(jj, 0) - 2 * (gamaj_tob.Point(0, 0) / betaj) * Vk.Point(jj, 0));
				}
				flag = true;
			}

		}
	}
}
//ֻ��A��Ԫ
void Matrix::householder()
{
	Matrix Vk(rows_num, 1);
	Matrix aj(rows_num, 1);
	for (int k = 0; k < cols_num; k++)
	{
		double alphak = -sgn(this->Point(k, k)) * (this->col_norm_2(k, k));
		//��Vk��ֵ
		for (int i = 0; i < rows_num; i++)
		{
			if (i < k)
			{
				Vk.setpoint(i, 0, 0);
			}
			else if (i == k)
			{
				Vk.setpoint(i, 0, (this->Point(k, k)) - alphak);
			}
			else if (i > k)
			{
				Vk.setpoint(i, 0, this->Point(i, k));
			}
		}
		double betaj = Vk.norm2(0);
		if (betaj == 0)
		{
			continue;
		}
		for (int j = k; j < cols_num; j++)
		{
			//����aj�ĸ���
			for (int ii = 0; ii < rows_num; ii++)
			{
				aj.setpoint(ii, 0, this->Point(ii, j));
			}
			Matrix gamaj(1, 1, 0);
			Matrix Vkt(Vk.col(), Vk.row());
			Vkt = Matrix::T(Vk);
			gamaj = Vkt * aj;
			//�õ�aj
			for (int jj = 0; jj < rows_num; jj++)
			{
				this->setpoint(jj, j, this->Point(jj, j) - 2 * (gamaj.Point(0, 0) / betaj) * Vk.Point(jj, 0));
			}
		}
	}
}
//ĳһ�е�k��Ԫ�غ��2����
double Matrix::col_norm_2(int current_col, int k)
{
	double sum = 0;
	for (int i = k; i < rows_num; i++)
	{
		sum += pow(this->Point(i, current_col), 2);
	}
	return sqrt(sum);
}
//���ź���
int sgn(double d) { return (d < 0 ? (-1) : (1)); }
//����Ԫ�ص�ֵ
void Matrix::setpoint(int row, int col, double m)
{
	this->p[row][col] = m;
}
//ĳһ�е�2������ƽ��
double Matrix::norm2(int which_col)
{
	double sum = 0;
	for (int i = 0; i < rows_num; i++)
	{
		sum += pow(this->Point(i, which_col), 2);
	}
	return sum;
}
//���ƹ��캯��
Matrix::Matrix(const Matrix& A)
{
	p = new double* [A.rows_num];
	for (int i = 0; i < A.rows_num; i++)
	{
		p[i] = new double[A.cols_num];
	}
	for (int i = 0; i < A.rows_num; i++)
	{
		for (int j = 0; j < A.cols_num; j++)
		{
			p[i][j] = A.Point(i, j);
		}
	}
	rows_num = A.rows_num;
	cols_num = A.cols_num;
}
//��ȡ�Ӿ���
Matrix Matrix::partial_Matrix(int row_begin, int row_end, int col_begin, int col_end)
{
	Matrix A(row_end - row_begin + 1, col_end - col_begin + 1);
	int ii = 0; int jj = 0;
	for (int i = row_begin - 1; i < row_end;)
	{
		for (int j = col_begin - 1; j < col_end; )
		{
			A.setpoint(ii, jj, this->Point(i, j));
			j++;
			if (j < col_end)
				jj++;
		}
		jj = 0;
		i++;
		if (i < row_end)
			ii++;
	}
	return A;
}
//��������
Matrix Matrix::number_times(double m)
{
	Matrix result(this->rows_num, this->cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			result.setpoint(i, j, this->Point(i, j) * m);
		}
	}
	return result;
}
//���淽�������
Matrix Matrix::normalization(Matrix& A, Matrix& b)
{
	Matrix x(A.cols_num, 1);
	Matrix At = Matrix::T(A);
	//x=Matrix::Solve(At*A,At*b);
	return x;
}
//���㻯
Matrix Matrix::augmented_solve(const Matrix& A, const Matrix& b)
{
	int m = A.rows_num;
	int n = A.cols_num;
	int a = m + n;
	Matrix augment(a, a);
	Matrix b_1(a, 1);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j)
			{
				augment.setpoint(i, j, 1);
			}//������λ��
			else
				augment.setpoint(i, j, 0);
		}
	}
	int ii = 0;
	int jj = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = m; j < a; j++)
		{
			augment.setpoint(i, j, A.Point(ii, jj));//��A����
			jj++;
		}
		jj = 0;
		ii++;
	}
	Matrix AT(n, m);
	AT = Matrix::T(A);
	ii = 0;
	jj = 0;
	for (int i = m; i < a; i++)
	{
		for (int j = 0; j < m; j++)
		{
			augment.setpoint(i, j, AT.Point(ii, jj));//��A��ת�þ���
			jj++;
		}
		jj = 0;
		ii++;
	}
	for (int i = m; i < a; i++)
	{
		for (int j = m; j < a; j++)
		{
			augment.setpoint(i, j, 0);//��A��ת�þ���
		}
	}
	for (int i = 0; i < a; i++)
	{
		if (i < m)
			b_1.setpoint(i, 0, b.Point(i, 0));
		else
			b_1.setpoint(i, 0, 0);
	}
	Matrix x(a, 1);
	x = Matrix::Solve(augment, b_1);
	return x;
}
void Matrix::givens_rotation(Matrix& A, Matrix& b)
{

	for (int i = 0; i < A.cols_num; i++)//�ӵ�һ����ʼ��Ԫ
	{
		int position = Matrix::find_nonzero(A, i);
		while (position != 0)
		{
			Matrix givens(A.rows_num, A.rows_num);
			for (int k = 0; k < A.rows_num; k++)//������λ��
			{
				givens.setpoint(k, k, 1);
			}
			double cos = A.Point(i, i) / sqrt(pow(A.Point(position, i), 2) + pow(A.Point(i, i), 2));
			double sin = A.Point(position, i) / sqrt(pow(A.Point(position, i), 2) + pow(A.Point(i, i), 2));
			givens.setpoint(i, i, cos);
			givens.setpoint(i, position, sin);
			givens.setpoint(position, i, -sin);
			givens.setpoint(position, position, cos);
			//�����Ԫ
			A = givens * A;
			b = givens * b;
			//A.Show();
			//b.Show();
			position = Matrix::find_nonzero(A, i);

		}
	}
}
int Matrix::find_nonzero(const Matrix& A, int which_col)//����ĳһ�еķ���Ԫ��
{
	bool ifreturn = false;
	for (int row = A.rows_num - 1; row > which_col; row--)//�������
	{
		if (abs(A.Point(row, which_col)) != 0)
		{
			return row;
			ifreturn = true;
			break;
		}
	}
	if (ifreturn == false)
		return 0;
}

Matrix* Matrix::gram_schmidt(Matrix& A, Matrix& b)//gram schmidt�任
{
	Matrix Q(A.rows_num, A.cols_num);
	Matrix R(A.cols_num, A.cols_num);
	/*Matrix qj(A.rows_num, 1);
	Matrix ak(A.rows_num, 1);
	for (int k = 0; k < A.cols_num; k++)
	{
		for (int j = 0; j <k - 1; j++)
		{
			qj = Q.partial_Matrix(1, A.rows_num, j+1, j+1);
			ak = A.partial_Matrix(1, A.rows_num, k+1, k+1);
			Matrix temp(1, 1);
			temp=Matrix::T(qj)* ak;
			R.setpoint(j, k, temp.Point(0, 0));
			for (int m = 0; m < A.rows_num; m++)
			{
				Q.setpoint(m, k, Q.Point(m, k) - temp.Point(0, 0) * qj.Point(m, 0));
			}	//jianche xiangguan
		}
		R.setpoint(k, k, sqrt(Q.norm2(k)));
		for (int m = 0; m < A.rows_num; m++)
		{
			Q.setpoint(m, k, Q.Point(m,k)/ sqrt(Q.norm2(k)));
		}
	}
	static Matrix result[2] = { Q,R };
	result[0] = Q;
	result[1] = R;
	return result;
	*/
	for (int k = 0; k < A.cols_num; k++)
	{
		Matrix qk(A.rows_num, 1);//qk=ak
		for (int i = 0; i < A.rows_num; i++)
		{
			qk.setpoint(i, 0, A.Point(i, k));
		}
		for (int j = 0; j <= k - 1; j++)
		{
			//�õ�qj:
			Matrix qj(A.rows_num, 1);//qk=ak
			for (int i = 0; i < A.rows_num; i++)
			{
				qj.setpoint(i, 0, Q.Point(i, j));
			}
			Matrix temp(1, 1);
			temp = Matrix::T(qj) * qk;
			double rjk = temp.Point(0, 0);
			R.setpoint(j, k, rjk);
			//qj=qj.number_times(rjk);//qj����rjk
			qk = (qk - (qj * rjk));
		}
		double rkk = sqrt(qk.norm2(0));
		R.setpoint(k, k, rkk);
		if (rkk == 0)
		{
			cout << " ����������Ϊ0��ֹͣ����\n";
			break;
		}
		else
		{
			qk = qk / rkk;
			for (int i = 0; i < A.rows_num; i++)//����ak
			{
				Q.setpoint(i, k, qk.Point(i, 0));
			}
		}
	}
	static Matrix result[2] = { Q,R };
	result[0] = Q;
	result[1] = R;
	return result;

}
Matrix* Matrix::gram_improve(Matrix& A, Matrix& b)
{
	Matrix Q(A.rows_num, A.cols_num);
	Matrix R(A.cols_num, A.cols_num);
	for (int k = 0; k < A.cols_num; k++)
	{
		double rkk = sqrt(A.norm2(k));
		R.setpoint(k, k, rkk);
		if (rkk == 0)
		{
			cout << " ����������Ϊ0��ֹͣ����\n";
			break;
		}
		Matrix qk(A.rows_num, 1);// jiang qk biao zhun hua
		for (int i = 0; i < A.rows_num; i++)
		{
			qk.setpoint(i, 0, A.Point(i, k) / rkk);//biao zhun hua
		}

		for (int i = 0; i < A.rows_num; i++)
		{
			Q.setpoint(i, k, qk.Point(i, 0));
		}
		for (int j = k + 1; j < A.cols_num; j++)
		{

			Matrix aj(A.rows_num, 1);
			for (int i = 0; i < A.rows_num; i++)
			{
				aj.setpoint(i, 0, A.Point(i, j));
			}
			Matrix temp(1, 1);
			temp = Matrix::T(qk) * aj;
			double rkj = temp.Point(0, 0);
			R.setpoint(k, j, rkj);
			aj -= qk.number_times(rkj);
			for (int i = 0; i < A.rows_num; i++)
			{
				A.setpoint(i, j, aj.Point(i, 0));
			}
		}
	}
	static Matrix result[2] = { Q,R };
	result[0] = Q;
	result[1] = R;
	return result;

}
Matrix* Matrix::iteration_gram_schmidt(Matrix& A, Matrix& b)//iteration
{
	//��һ�ηֽ���ľ���
	Matrix tempQ(A.rows_num, A.cols_num);
	Matrix tempR(A.cols_num, A.cols_num);
	Matrix* temp;
	temp = Matrix::gram_schmidt(A, b);
	tempQ = temp[0];
	tempR = temp[1];
	Matrix* answer = Matrix::gram_improve(tempQ, b);
	answer[1] = answer[1] * tempR;
	return answer;
}
Matrix Matrix::operator/(double m)
{
	Matrix re(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			re.setpoint(i, j, this->Point(i, j) / m);
		}
	}
	return re;
}
double Matrix::infinite_norm(int which_col)
{
	double max = 0;
	for (int i = 0; i < rows_num; i++)
	{
		if (abs(Point(i, which_col)) > max)
		{
			max = abs(Point(i, which_col));
		}
	}
	return max;
}
Matrix Matrix::operator+(const Matrix& m)
{
	Matrix temp(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			temp.p[i][j] = p[i][j] + m.p[i][j];
		}
	}
	return temp;
}

Matrix Matrix::operator-(Matrix m)
{
	Matrix temp(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			temp.p[i][j] = p[i][j] - m.p[i][j];
		}
	}
	return temp;
}
Matrix Matrix::operator*(double m)
{
	Matrix re(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			re.setpoint(i, j, this->Point(i, j) * m);
		}
	}
	return re;
}
void Matrix::gaussianEliminate(Matrix& b)
{
	Matrix Ab(*this);
	int rows = Ab.rows_num;
	int cols = Ab.cols_num;
	int Acols = cols - 1;

	int i = 0; //������
	int j = 0; //������
	while (i < rows)
	{
		bool flag = false;
		while (j < Acols && !flag)
		{
			if (Ab.p[i][j] != 0)
			{
				flag = true;
			}
			else
			{
				int max_row = i;
				double max_val = 0;
				for (int k = i + 1; k < rows; ++k)
				{
					double cur_abs = Ab.p[k][j] >= 0 ? Ab.p[k][j] : -1 * Ab.p[k][j];
					if (cur_abs > max_val)
					{
						max_row = k;
						max_val = cur_abs;
					}
				}
				if (max_row != i)
				{
					Ab.swapRows(max_row, i);
					b.swapRows(max_row, i);
					flag = true;
				}
				else
				{
					j++;
				}
			}
		}
		if (flag)
		{
			for (int t = i + 1; t < rows; t++)
			{
				for (int s = j + 1; s < cols; s++)
				{
					Ab.p[t][s] = Ab.p[t][s] - Ab.p[i][s] * (Ab.p[t][j] / Ab.p[i][j]);

					if (abs(Ab.p[t][s]) <= EPS)
						Ab.p[t][s] = 0;
				}
				Ab.p[t][0] = b.p[t][0] - b.p[i][0] * (Ab.p[t][j] / Ab.p[i][j]);
				Ab.p[t][j] = 0;
			}
		}
		i++;
		j++;
	}
	*this = Ab;
}
Matrix Matrix::Solve_only(Matrix& A, Matrix& b)
{
	Matrix x(b.rows_num, 1);
	//�������һ��
	x.p[x.rows_num - 1][0] = b.p[x.rows_num - 1][0] / A.p[x.rows_num - 1][x.rows_num - 1];

	if (abs(x.p[x.rows_num - 1][0]) < EPS)
		x.p[x.rows_num - 1][0] = 0;
	for (int i = x.rows_num - 2; i >= 0; i--)//���лش�
	{
		double sum = 0;
		for (int j = i + 1; j < x.rows_num; j++)
		{
			sum += A.p[i][j] * x.p[j][0];
		}
		x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
		if (abs(x.p[i][0]) < EPS)
			x.p[i][0] = 0;
	}

	return x;
}
int Matrix::find_max(const Matrix& A, int which_col)
{
	double max = 0;
	int re = 0;
	for (int row = A.rows_num - 1; row >= which_col; row--)//�������
	{
		if (abs(A.Point(row, which_col)) > max)
		{
			max = abs(A.Point(row, which_col));
			re = row;
		}
	}
	return re;
}
Matrix* Matrix::gram_schmidt(Matrix& A)//gram�����任
{
	Matrix Q(A.rows_num, A.cols_num);
	Matrix R(A.cols_num, A.cols_num);
	/*Matrix q1(A.rows_num,A.cols_num);
	Matrix q[10];

	for (int i = 0; i < A.rows_num; i++)
	{
		q[i].p=new double* [A.rows_num];//����rows_num��ָ��
		for (int m = 0; m < A.rows_num; ++m)
		{
			q[i].p[m] = new double[A.cols_num];//Ϊp[i]���ж�̬�ڴ���䣬��СΪcols
		}
		q[i].rows_num = A.rows_num;
		q[i].cols_num = A.cols_num;
		q[i] = Matrix(A.rows_num, 1);
	}
	for (int k = 0; k < A.cols_num; k++)
	{
		for (int i = 0; i < A.rows_num; i++)
		{
			q[k].setpoint(i, 0, A.Point(i, k));
		}
		for (int j = 0; j <= k - 1; j++)  // �ӵ�ǰ���м�ȥ��ǰ���еķ���
		{
			R.setpoint(j, k, (Matrix::T(q[j]) * q[k]).Point(0, 0));
			//R.MT[j][k] = (!q[j] * q[k]).MT[0][0];

			for (int i = 0; i < A.rows_num; i++)
			{
				q[k].setpoint(i, 0, q[k].Point(i, 0) - R.Point(j, k) * q[j].Point(i, 0));
				//q[k].MT[i][0] = q[k].MT[i][0] - R.MT[j][k] * q[j].MT[i][0];
			}
		}
		R.setpoint(k, k, sqrt((Matrix::T(q[k]) * q[k]).Point(0, 0)));
		for (int i = 0; i < A.rows_num; i++)	// ����ǰ�б�׼��
		{
			Q.p[i][k] = q[k].p[i][0] = q[k].p[i][0] / R.p[k][k];
		}

	}
	*/

	for (int k = 0; k < A.cols_num; k++)
	{
		Matrix qk(A.rows_num, 1);//qk=ak
		for (int i = 0; i < A.rows_num; i++)
		{
			qk.setpoint(i, 0, A.Point(i, k));
		}
		cout << "qk" << endl;
		qk.Show();
		for (int j = 0; j <= k - 1; j++)
		{
			//�õ�qj:
			Matrix qj(A.rows_num, 1);
			for (int i = 0; i < A.rows_num; i++)
			{
				qj.setpoint(i, 0, Q.Point(i, j));
			}
			cout << "qj" << endl;
			qj.Show();
			Matrix temp(1, 1);
			temp = Matrix::T(qj) * qk;
			double rjk = temp.Point(0, 0);
			cout << "rjk" << endl;
			cout << rjk << endl;
			R.setpoint(j, k, rjk);
			//qj = qj.number_times(rjk);//qj����rjk
			qk -= qj * rjk;
			cout << "qk" << endl;
			qk.Show();
		}
		double rkk = sqrt(qk.norm2(0));
		cout << "rkk" << endl;
		cout << rkk << endl;
		R.setpoint(k, k, rkk);
		if (rkk == 0)
		{
			cout << " ����������Ϊ0��ֹͣ����\n";
			break;
		}
		else
		{
			qk = qk / rkk;
			cout << "qk" << endl;
			qk.Show();
			for (int i = 0; i < A.rows_num; i++)//����ak
			{
				Q.setpoint(i, k, qk.Point(i, 0));
			}
		}
	}

	static Matrix result[2] = { Q,R };
	result[0] = Q;
	result[1] = R;
	//Q.Show();
	//R.Show();
	return result;
}
long long factorial(long long n)//�׳˺���
{
	if (n == 1)
	{
		return 1;
	}
	else
	{
		return n * factorial(n - 1);
	}
}
Matrix Matrix::exp_sum(int n)
{
	Matrix re(rows_num, cols_num);
	if (n == 0)
	{
		return Matrix::eye(rows_num);
	}
	else
	{
		re = (((*this) ^ n) / factorial(n)) + this->exp_sum(n - 1);
	}
	return re;
}

Matrix Matrix::operator^(int n)//����^
{
	Matrix temp(rows_num, cols_num);
	if (n == 1)
	{
		temp = *this;
	}
	else
	{
		temp = (*this) * (*this ^ (n - 1));
	}
	return temp;
}
Matrix Matrix::exp_gigenvector()
{
	//Matrix* QR;
	//Matrix b(rows_num, 1);
	//Matrix pre = *this();

	Matrix Q1(rows_num, cols_num);
	for (int i = 0; i < 10; i++)
	{

		Matrix Q(rows_num, cols_num);
		Matrix R(cols_num, cols_num);
		Matrix::gram_schmidt(*this, Q, R);
		*this = (R * Q);
		//cout << "qr" << endl;
		//QR[0].Show();
		//QR[1].Show();
		//cout << "xiang cheng" << endl;
		//(QR[0]*QR[1]).Show();
		//cout << "this" << endl;
		//this->Show();

		if (i == 0)
		{
			Q1 = Q;
		}
		else
		{
			Q1 = (Q1 * Q);
		}
	}

	//Q��������ΪA����������
	Matrix D(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		D.setpoint(i, i, exp(this->Point(i, i)));
	}
	Q1.Show();
	Matrix re(rows_num, cols_num);
	re = Q1 * D * (Matrix::inv(Q1));
	return re;
}
Matrix* Matrix::gram_improve(Matrix& A)
{
	Matrix Q(A.rows_num, A.cols_num);
	Matrix R(A.cols_num, A.cols_num);
	for (int k = 0; k < A.cols_num; k++)
	{
		double rkk = sqrt(A.norm2(k));
		//cout << rkk << endl;
		R.setpoint(k, k, rkk);
		if (rkk == 0)
		{
			cout << " ����������Ϊ0��ֹͣ����\n";
			break;
		}
		Matrix qk(A.rows_num, 1);// jiang qk biao zhun hua
		for (int i = 0; i < A.rows_num; i++)
		{
			qk.setpoint(i, 0, A.Point(i, k) / rkk);//biao zhun hua
		}
		//cout << "qk" << endl;
		//qk.Show();
		for (int i = 0; i < A.rows_num; i++)
		{
			Q.setpoint(i, k, qk.Point(i, 0));
		}
		for (int j = k + 1; j < A.cols_num; j++)
		{

			Matrix aj(A.rows_num, 1);
			for (int i = 0; i < A.rows_num; i++)
			{
				aj.setpoint(i, 0, A.Point(i, j));
			}
			Matrix temp(1, 1);
			temp = Matrix::T(qk) * aj;
			double rkj = temp.Point(0, 0);
			//	cout << "rkj" << endl;
			//	cout << rkj << endl;
			R.setpoint(k, j, rkj);
			aj = aj - qk * rkj;
			for (int i = 0; i < A.rows_num; i++)
			{
				A.setpoint(i, j, aj.Point(i, 0));
			}
			//for (int i = 0; i < A.rows_num; i++)
			//{
				//A.setpoint(i, j, qk.Point(i, 0));
		///	}
		//	cout << "A" << endl;
		//	A.Show();
		}
	}
	static Matrix result[2] = { Q,R };
	result[0] = Q;
	result[1] = R;
	return result;

}
void Matrix::gram_schmidt(Matrix& A, Matrix& Q, Matrix& R)//gram�����任
{
	//Matrix Q(A.rows_num, A.cols_num);
	//Matrix R(A.cols_num, A.cols_num);
	//Matrix q1(A.rows_num,A.cols_num);
	Matrix q[10];

	for (int i = 0; i < A.cols_num; i++)
	{
		q[i] = Matrix(A.rows_num, 1);
	}
	for (int k = 0; k < A.cols_num; k++)
	{
		for (int i = 0; i < A.rows_num; i++)
		{
			q[k].setpoint(i, 0, A.Point(i, k));
		}
		for (int j = 0; j <= k - 1; j++)  // �ӵ�ǰ���м�ȥ��ǰ���еķ���
		{
			R.setpoint(j, k, (Matrix::T(q[j]) * q[k]).Point(0, 0));
			//R.MT[j][k] = (!q[j] * q[k]).MT[0][0];

			for (int i = 0; i < A.rows_num; i++)
			{
				q[k].setpoint(i, 0, q[k].Point(i, 0) - R.Point(j, k) * q[j].Point(i, 0));
				//q[k].MT[i][0] = q[k].MT[i][0] - R.MT[j][k] * q[j].MT[i][0];
			}
		}
		R.setpoint(k, k, sqrt((Matrix::T(q[k]) * q[k]).Point(0, 0)));
		for (int i = 0; i < A.rows_num; i++)	// ����ǰ�б�׼��
		{
			Q.p[i][k] = q[k].p[i][0] = (q[k].p[i][0] / R.p[k][k]);
		}

	}

}
Matrix Matrix::diag()//shuangcexuanzhuan
{
	Matrix a = *this;
	if (this->Point(0, 1) != this->Point(1, 0))
	{
		cout << "���Գ�����" << endl;

	}
	else
	{
		Matrix Q1(rows_num, cols_num);
		for (int i = 0; i < 10; i++)
		{

			Matrix Q(rows_num, cols_num);
			Matrix R(cols_num, cols_num);
			Matrix::gram_schmidt(*this, Q, R);
			*this = (R * Q);

			if (i == 0)
			{
				Q1 = Q;
			}
			else
			{
				Q1 = (Q1 * Q);
			}

		}
		*this = a;
		return Q1;
	}
}
Matrix Matrix::res(int which_row, int which_col)
{
	Matrix temp(rows_num - 1, cols_num - 1);
	int k = 0;
	for (int temp_row = 0; temp_row < temp.rows_num; temp_row++)
	{
		if (temp_row < which_row)
		{
			for (int temp_col = 0; temp_col < temp.cols_num; temp_col++)
			{
				if (temp_col < which_col)
				{
					temp.setpoint(temp_row, temp_col, this->Point(k, temp_col));
				}
				else
				{
					temp.setpoint(temp_row, temp_col, this->Point(k, temp_col + 1));
				}
			}
		}
		else
		{
			for (int temp_col = 0; temp_col < temp.cols_num; temp_col++)
			{
				if (temp_col < which_col)
				{
					temp.setpoint(temp_row, temp_col, this->Point(k + 1, temp_col));
				}
				else
				{
					temp.setpoint(temp_row, temp_col, this->Point(k + 1, temp_col + 1));
				}
			}
		}
		k++;
		
	}
	return temp;
}
Matrix Matrix::combine_rows(Matrix& A,Matrix& B)
{
	if (A.cols_num == B.cols_num)
	{
		Matrix result(A.rows_num + B.rows_num, A.cols_num);
		for (int j = 0; j < A.cols_num; j++)
		{
			for (int i = 0; i < A.rows_num; i++)
			{
				result.p[i][j] = A.p[i][j];
			}
			for (int i = A.rows_num; i < result.rows_num; i++)
			{
				result.p[i][j] = B.p[i - A.rows_num][j];
			}
		}
		return result;
	}
	else
	{
		cout << "�������ȣ������кϲ�" << endl;
		abort();
	}
}
Matrix Matrix::combine_cols(Matrix& A, Matrix& B)
{
	if (A.rows_num == B.rows_num)
	{
		Matrix result(A.rows_num, A.cols_num + B.cols_num);
		for (int i = 0; i < A.rows_num; i++)
		{
			for (int j = 0; j < A.cols_num; j++)
			{
				result.p[i][j] = A.p[i][j];
			}
			for (int j = A.cols_num; j < result.cols_num; j++)

			{
				result.p[i][j] = B.p[i][j - A.cols_num];
			}
		}
		return result;
	}
	else
	{
		cout << "�������ȣ������кϲ�" << endl;
		abort();
	}
}
