#include "complex_matrix.h"
#include"complex.h"
#include<iomanip>
void complex_matrix::initialize()
{
	p = new Complex *[rows_num];//分配rows_num个指针
	for (int i = 0; i < rows_num; ++i)
	{
		p[i] = new Complex[cols_num];//为p[i]进行动态内存分配，大小为cols
	}
}
complex_matrix::complex_matrix(int row, int col)
{
	rows_num = row;
	cols_num = col;
	initialize();
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			p[i][j].real =0;
			p[i][j].imag =0;//初始化
		}
	}
}
complex_matrix::~complex_matrix() 
{
	for (int i = 0; i < rows_num; ++i)
	{
		delete[] p[i];
	}
	delete[] p;
}
complex_matrix::complex_matrix()
{
	rows_num = 1;
	cols_num = 1;
	initialize();
}
complex_matrix::complex_matrix(const complex_matrix&A)
{
	p = new Complex* [A.rows_num];
	for (int i = 0; i < A.rows_num; i++)
	{
		p[i] = new Complex [A.cols_num] ;
	}
	for (int i = 0; i < A.rows_num; i++)
	{
		for (int j = 0; j < A.cols_num; j++)
		{
			p[i][j] = A.p[i][j];
		}
	}
	rows_num = A.rows_num;
	cols_num = A.cols_num;
}
complex_matrix complex_matrix::operator=(const complex_matrix& m)
{
	if (this == &m)
	{
		return *this;
	}
	if (rows_num != m.rows_num || cols_num != m.cols_num)
	{
		//cout << "维度不相同" << endl;
		for (int i = 0; i < rows_num; ++i) 
		{
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
Complex complex_matrix::operator()(int row, int col)const//返回在某一点的值
{
	return this->p[row][col];
}
void complex_matrix::Show() const//矩阵显示
{
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			cout << setw(7) << setprecision(6) << fixed << " " << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
complex_matrix& complex_matrix::operator=(double* a)
{
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			p[i][j].real =* (a + i * cols_num + j);
			p[i][j].imag = 0;
		}
	}
	return *this;
}
complex_matrix& complex_matrix::operator=(Complex* a)
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
complex_matrix complex_matrix::operator*(const complex_matrix& m)
{
	if (cols_num != m.rows_num)
	{
		cout << "维度不等" << endl;
	}
	complex_matrix ba_M(rows_num, m.cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < m.cols_num; j++)
		{
			for (int k = 0; k < cols_num; k++)
			{
				Complex temp;
				temp = (p[i][k] * m.p[k][j]);
				ba_M.p[i][j] = ba_M.p[i][j] + (p[i][k] * m.p[k][j]);
			}
		}
	}
	return ba_M;
}
complex_matrix complex_matrix::eyes(int n)
{
	complex_matrix result(n, n);
	for (int i = 0; i < n; i++)
	{
		result.p[i][i].real = 1.0;
	}
	return result;
}
complex_matrix complex_matrix::operator+(const complex_matrix& A)
{
	if (this->rows_num != A.rows_num || this->cols_num != A.cols_num)	//维度不匹配则不能运算
	{
		cout << "维度不匹配，不能执行矩阵加法" << endl;
		system("pause");
	}
	complex_matrix re(this->rows_num, this->cols_num);
	for (int i = 0; i < re.rows_num; i++)
	{
		for (int j = 0; j < re.cols_num; j++)
		{
			re.p[i][j] = this->p[i][j] + A.p[i][j];
		}
	}
	return re;
}
complex_matrix complex_matrix::operator-(const complex_matrix& A)
{
	if (this->rows_num != A.rows_num || this->cols_num != A.cols_num)	//维度不匹配则不能运算
	{
		cout << "维度不匹配，不能执行矩阵加法" << endl;
		system("pause");
	}
	complex_matrix re(this->rows_num, this->cols_num);
	for (int i = 0; i < re.rows_num; i++)
	{
		for (int j = 0; j < re.cols_num; j++)
		{
			re.p[i][j] = this->p[i][j] - A.p[i][j];
		}
	}
	return re;
}
complex_matrix complex_matrix::operator*(double m)
{
	complex_matrix re(this->rows_num, this->cols_num);
	for (int i = 0; i < re.rows_num; i++)
	{
		for (int j = 0; j < re.cols_num; j++)
		{
			re.p[i][j] = this->p[i][j] *m;
		}
	}
	return re;
}
complex_matrix complex_matrix::operator/(double m)
{
	complex_matrix re(this->rows_num, this->cols_num);
	for (int i = 0; i < re.rows_num; i++)
	{
		for (int j = 0; j < re.cols_num; j++)
		{
			re.p[i][j] = this->p[i][j] / m;
		}
	}
	return re;
}
complex_matrix complex_matrix::operator!()const
{
	complex_matrix re(this->cols_num, this->rows_num);
	for (int i = 0; i < re.rows_num; i++)
	{
		for (int j = 0; j < re.cols_num; j++)
		{
			re.p[i][j] = !this->p[j][i];
		}
	}
	return re;
}
double complex_matrix::norm2(int which_col)
{
	double sum = 0;
	for (int i = 0; i < rows_num; i++)
	{
		sum += pow(p[i][which_col].modulus(), 2);
	}
	return sqrt(sum);
}
complex_matrix complex_matrix::col_vector(int which_col)
{
	complex_matrix temp(rows_num, 1);
	for (int i = 0; i < rows_num; i++)
	{
		temp.p[i][0] = this->p[i][which_col];
	}
	return temp;
}
complex_matrix complex_matrix::row_vector(int which_row)
{
	complex_matrix temp(1, cols_num);
	for (int i = 0; i < cols_num; i++)
	{
		temp.p[0][i] = this->p[which_row][i];
	}
	return temp;
}
int complex_matrix::row() const
{
	return rows_num;
}
int complex_matrix::col() const
{
	return cols_num;
}

complex_matrix complex_matrix::inverse()//使用拆分方式求逆,方阵
{
	Matrix A = this->realmatrix();///实部矩阵
	Matrix B = this->imagmatrix();
	Matrix tempreal(rows_num, cols_num);//逆矩阵的实部矩阵
	Matrix tempimag(rows_num, cols_num);
	complex_matrix inv(rows_num, cols_num);//逆矩阵
	tempreal = Matrix::inv(A + B * Matrix::inv(A) * B);
	tempimag = Matrix::inv(A) * B * (tempreal);
	tempimag = Matrix(rows_num, cols_num, 0.0) - tempimag;
	inv.set_real_matrix(tempreal);
	inv.set_imag_matrix(tempimag);
	return inv;
}
Matrix complex_matrix::realmatrix()
{
	Matrix temp(this->rows_num, this->cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			temp.setpoint(i, j, this->p[i][j].real);
		}
	}
	return temp;
}
Matrix complex_matrix::imagmatrix()
{
	Matrix temp(this->rows_num, this->cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			temp.setpoint(i, j, this->p[i][j].imag);
		}
	}
	return temp;
}
complex_matrix& complex_matrix::set_real_matrix(const Matrix& A)
{
	//complex_matrix temp(this->rows_num, this->cols_num);
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			this->p[i][j].real = A.Point(i, j);
		}
	}
	return *this;
}
complex_matrix& complex_matrix::set_imag_matrix(const Matrix& A)
{
	
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
		{
			this->p[i][j].imag = A.Point(i, j);
		}
	}
	return *this;
}
complex_matrix complex_matrix::pseudo_inverse()const//伪逆
{
	complex_matrix A=*this;
	complex_matrix A_H = !A;
	complex_matrix temp = A_H * A;
	temp = temp.inverse();
	temp = temp * A_H;
	return temp;
}
//得到第i行
complex_matrix complex_matrix::get_row(int i)
{
	complex_matrix target_row(1, cols_num);
	for (int k = 0; k < cols_num; k++) target_row.p[0][k] = p[i][k];
	return target_row;
}
//得到第j列
complex_matrix complex_matrix::get_col(int j)
{
	complex_matrix target_column(rows_num,1);
	for (int k = 0; k < rows_num; k++) target_column.p[k][0] = p[k][j];
	return target_column;
}
double complex_matrix::norm2_vector()//向量的2范数
{
	double result = 0.0;
	if (rows_num == 1)		//行向量
	{
		for (int j = 0; j < cols_num; j++)
			result = result + p[0][j].real * p[0][j].real + p[0][j].imag * p[0][j].imag;
		result = sqrt(result);
		return result;
	}
	else if (cols_num == 1)	//列向量
	{
		for (int i = 0; i < rows_num; i++)
			result = result + p[i][0].real * p[i][0].real + p[i][0].imag * p[i][0].imag;
		result = sqrt(result);
		return result;
	}
}
complex_matrix complex_matrix::remove_row(int which_row)
{
	complex_matrix temp(this->rows_num - 1, this->cols_num);
	for (int row = 0; row < temp.rows_num; row++)
	{
		for (int col = 0; col < temp.cols_num; col++)
		{
			if (row < which_row)
			{
				temp.p[row][col] = this->p[row][col];
			}
			else
			{
				temp.p[row][col] = this->p[row + 1][col];
			}
		}
	}
	return temp;
}
complex_matrix complex_matrix::remove_col(int which_col)
{
	complex_matrix temp(this->rows_num , this->cols_num-1);
	for (int row = 0; row < temp.rows_num; row++)
	{
		for (int col = 0; col < temp.cols_num; col++)
		{
			if (col < which_col)
			{
				temp.p[row][col] = this->p[row][col];
			}
			else
			{
				temp.p[row][col] = this->p[row ][col+1];
			}
		}
	}
	return temp;
}
complex_matrix complex_matrix::operator*(Complex c)
{
	complex_matrix temp = *this;
	for (int i = 0; i < rows_num; i++)
	{
		for (int j = 0; j < cols_num; j++)
			temp.p[i][j] = temp.p[i][j] * c;
	}
	return temp;
}
complex_matrix* complex_matrix::gram_improve(complex_matrix& A)
{
	complex_matrix temp = A;
	complex_matrix Q(A.rows_num, A.cols_num);
	complex_matrix R(A.cols_num, A.cols_num);
	for (int k = 0; k < A.cols_num; k++)
	{
		double rkk = A.get_col(k).norm2_vector();
		R.p[k][k].real=rkk;

		if (rkk == 0)
		{
			cout << " 矩阵列向量为0，停止计算\n";
			break;
		}
		complex_matrix qk(A.rows_num, 1);// jiang qk biao zhun hua
		for (int i = 0; i < A.rows_num; i++)
		{
			qk.p[i][0] = (A.p[i][k] / rkk);// setpoint(i, 0, A.Point(i, k) / rkk);//biao zhun hua
		}
		//cout << "qk" << endl;
		//qk.Show();
		for (int i = 0; i < A.rows_num; i++)
		{
			Q.p[i][k] = qk.p[i][0];// setpoint(i, k, qk.Point(i, 0));
		}
		for (int j = k + 1; j < A.cols_num; j++)
		{

			complex_matrix aj(A.rows_num, 1);
			for (int i = 0; i < A.rows_num; i++)
			{
				aj.p[i][0] = A.p[i][j];// setpoint(i, 0, A.Point(i, j));
			}
			complex_matrix temp(1, 1);
			temp = (!qk) * aj;
			Complex rkj = temp.p[0][0];// Point(0, 0);
			//	cout << "rkj" << endl;
			//	cout << rkj << endl;
			R.p[k][j] = rkj;// [setpoint(k, j, rkj);
			aj = aj - qk * rkj;
			for (int i = 0; i < A.rows_num; i++)
			{
				A.p[i][j] = aj.p[i][0];// setpoint(i, j, aj.Point(i, 0));
			}
			//for (int i = 0; i < A.rows_num; i++)
			//{
				//A.setpoint(i, j, qk.Point(i, 0));
		///	}
		//	cout << "A" << endl;
		//	A.Show();
		}
	}
	A = temp;
	static complex_matrix result[2] = { Q,R };
	result[0] = Q;
	result[1] = R;
	return result;
}
//换行
void complex_matrix::exchange_row(int i1, int i2)
{
	Complex temp;
	for (int j = 0; j < cols_num; j++)
	{
		temp = p[i1][j];
		p[i1][j] = p[i2][j];
		p[i2][j] = temp;
	}
}
void complex_matrix::exchange_col(int j1, int j2)
{	
	Complex temp;
	for (int i = 0; i < rows_num; i++)
	{
		temp = p[i][j1];
		p[i][j1] = p[i][j2];
		p[i][j2] = temp;
	}
}
complex_matrix complex_matrix::combine_rows(complex_matrix& A, complex_matrix& B)
{
	if (A.cols_num == B.cols_num)
	{
		complex_matrix result(A.rows_num + B.rows_num, A.cols_num);
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
		cout << "列数不等，不能行合并" << endl;
		abort();
	}
}
complex_matrix complex_matrix::combine_cols(complex_matrix& A, complex_matrix& B)
{
	if (A.rows_num == B.rows_num)
	{
		complex_matrix result(A.rows_num, A.cols_num + B.cols_num);
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
		cout << "行数不等，不能列合并" << endl;
		abort();
	}
}
Complex complex_matrix::cdot()
{
	if (cols_num == 1)//列向量情形
		return ((!(*this)) * (*this)).p[0][0];
	else if (rows_num == 1)
		return ((*this) * (!*this)).p[0][0];
	else
		cout << "不为向量请检查" << endl;
}
Complex complex_matrix::cdot(complex_matrix A)
{
	if (cols_num == 1)//列向量情形
		return ((!(*this)) * (A)).p[0][0];
	else if (rows_num == 1)
		return ((A) * (!*this)).p[0][0];
	else
		cout << "不为向量请检查" << endl;
}