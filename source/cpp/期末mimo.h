#pragma once
#include"matrix.h"
#include"complex_matrix.h"
#include"auxiliary.h"
#include<deque>
double f_s(Matrix s,Matrix A,Matrix b)//文献一中所用函数
{
	if (s.rows_num == A.rows_num && b.rows_num == s.rows_num)
		return ((Matrix::T(s)) * A * s * 0.5 - Matrix::T(b) * s).p[0][0];
	else
		cout << "维度不匹配请检查" << endl;

}
Complex f_s_BB(complex_matrix s, complex_matrix A, complex_matrix b,complex_matrix y)//bb算法中所用的格式
{
	if (s.rows_num == A.rows_num && b.rows_num == s.rows_num)
		return ((!s) * A * s * 0.5 - (!b) * s+(!y)*y*0.5).p[0][0];
	else
		cout << "维度不匹配请检查" << endl;
}




Matrix BFGS(complex_matrix H_hat,complex_matrix y_hat,double s,int iteration)//BFGS算法 s为方差 it 迭代次数
{
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y=R_y.combine_rows(R_y,I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H= H.combine_cols(H, temp);//最终得到H

	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;

	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix B = Matrix::eye(s0.rows_num);//初始时可以为任意对称正定矩阵
	Matrix g=get_negetive_gradient(*f_s,s0,A,b);//notice that g is negative 
	
	Matrix d;
	d = B * g;

	Matrix p;
	Matrix q;
	Matrix temps;
	Matrix tempg;
	
	double alpha = (Matrix::T(g) * d).p[0][0] / ((Matrix::T(d) * A * d)).p[0][0];
	//因为求得的g为负梯度，因此不需要加负号
	for (size_t i = 0; i < iteration; i++)
	{
		temps = s0;//save s0

		s0 = s0 + d*alpha;//迭代

		tempg = g;
		g = get_negetive_gradient(*f_s, s0, A, b);//update g, notice that g is negative

		//更新B
		p = s0 - temps;
		q = tempg-g;//because g is negative

		//q = Matrix(q.rows_num, q.cols_num) - q;//because g is negative. make it positive
		double ro = 1.0 / ((Matrix::T(p) * q).p[0][0]);
		Matrix V = Matrix::eye(s0.rows_num) - q * Matrix::T(p) * ro;
		B = Matrix::T(V) * B * V + p * Matrix::T(p) * ro;

		d = B * g;//更新dk

		alpha= (Matrix::T(g) * d).p[0][0] / ((Matrix::T(d) * A * d)).p[0][0];//更新α
			
	}
	return s0;
}
Matrix direct_inv(complex_matrix H_hat, complex_matrix y_hat,double s)//直接求逆的算法
{
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H

	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;
	Matrix re = Matrix::inv(A) * b;
	return re;
}
Matrix LBFGS_two_loop(int m, int k, deque<Matrix> deque_of_p
	, deque<Matrix> deque_of_q, deque<double> deque_of_ro, Matrix gk_1, Matrix Bk)//LBFBGS算法中的双向循环函数
{
	int delta;
	int L;
	if (k <= m)
	{
		 //delta = 0;
		 L = k;
	}
	else
	{
		//delta = k - m;
		L = m;
	}
	double* alpha = new double[L];
	Matrix r = gk_1;
	//int j;
	for (int i = L-1; i >=0; i--)
	{
		//j = i + delta;
		alpha[i] = (Matrix::T(deque_of_p[i]) * r * deque_of_ro[i]).p[0][0];
		r = r - deque_of_q[i]* alpha[i];
	} 
	r = Bk * r;
	double beta;
	for (int i = 0; i < L-1; i++)
	{
		beta = (Matrix::T( deque_of_q[i]) * r * deque_of_ro[i]).p[0][0];
		r = r + deque_of_p[i] * (alpha[i] - beta);
	}
	return Matrix(r.rows_num, r.cols_num) - r;
}
Matrix LBFGS(complex_matrix H_hat, complex_matrix y_hat, double s,int iteration,int m=5)//m默认为5 LBFGS
{
	deque<Matrix> deque_of_p;
	deque<Matrix> deque_of_q;
	deque<double> deque_of_ro;
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	//H.Show();
	//y.Show();


	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;



	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix g0 = A * s0 - b;
	Matrix B0 = Matrix::eye(s0.rows_num);//初始时可以为任意对称正定矩阵
	Matrix B;
	Matrix d0 = B0 * g0;
	d0 = Matrix(d0.rows_num, d0.cols_num) - d0;
	Matrix temps;
	Matrix tempg;
	Matrix p;
	Matrix q;
	for (int k = 0; k < iteration; k++)
	{
		/*if (g0.norm2(0) <= eps)
		{
			break;
		}*/
		//else
		{
			double alpha = -(Matrix::T(g0) * d0).p[0][0] / ((Matrix::T(d0) * A * d0)).p[0][0];

			temps = s0;//save s0

			s0 = s0 + d0 * alpha;//迭代

			tempg = g0;
			g0 = A * s0 - b;

			//更新B
			p = s0 - temps;
			q = g0 - tempg;
			double gamma = (Matrix::T(q) * p).p[0][0] / (Matrix::T(q) * q).p[0][0];
			B = B0 * gamma;
			//使用双端队列来保存pi和qi
			deque_of_p.push_back(p);
			deque_of_q.push_back(q);
			deque_of_ro.push_back(1.0 / (Matrix::T(p) * q).p[0][0]);
			if (k > m)
			{
				//k>m删除前端元素
				deque_of_p.pop_front();
				deque_of_q.pop_front();
				deque_of_ro.pop_front();
			}
			d0 = LBFGS_two_loop(m, k, deque_of_p, deque_of_q
				, deque_of_ro, g0, B);
		}
	}
	return s0;
}
Matrix LBFGS_B(complex_matrix H_hat, complex_matrix y_hat, double s, int iteration, int m = 5) //LBFGS算法中改进初始矩阵的策略
{
	deque<Matrix> deque_of_p;
	deque<Matrix> deque_of_q;
	deque<double> deque_of_ro;
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	//H.Show();
	//y.Show();


	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;

	

	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix g0 = A * s0 - b;
	Matrix B0 = Matrix::eye(s0.rows_num);//初始时可以为任意对称正定矩阵
	Matrix B;
	Matrix d0 = B0 * g0;
	d0 = Matrix(d0.rows_num, d0.cols_num) - d0;
	Matrix temps;
	Matrix tempg;
	Matrix p;
	Matrix q;
	for (int k = 0; k < iteration; k++)
	{
		/*if (g0.norm2(0)<=eps)
		{ 
			break;
		}
		else*/
		{
			double alpha = -(Matrix::T(g0) * d0).p[0][0] / ((Matrix::T(d0) * A * d0)).p[0][0];

			temps = s0;//save s0

			s0 = s0 + d0 * alpha;//迭代

			tempg = g0;
			g0 = A * s0 - b;

			//更新B
			p = s0 - temps;
			q = g0-tempg;
			double gamma = (Matrix::T(q)*p).p[0][0] / (Matrix::T(q)*q).p[0][0];
			B = B0 * gamma;
			//使用双端队列来保存pi和qi
			deque_of_p.push_back(p);
			deque_of_q.push_back(q);
			deque_of_ro.push_back(1.0 / (Matrix::T(p) * q).p[0][0]);
			if (k > m)
			{
				//k>m删除前端元素
				deque_of_p.pop_front();
				deque_of_q.pop_front();
				deque_of_ro.pop_front();
			}
			 d0 = LBFGS_two_loop(m,k,deque_of_p,deque_of_q
			,deque_of_ro,g0,B);
		}
	}
	return s0;
}
Matrix LBFGS_D(complex_matrix H_hat, complex_matrix y_hat, double s, int iteration, int m = 1)//LBFGS算法中 修改初始矩阵策略的一个分支
{
	deque<Matrix> deque_of_p;
	deque<Matrix> deque_of_q;
	deque<double> deque_of_ro;
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	//H.Show();
	//y.Show();


	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;
	Matrix D(A.rows_num, A.rows_num);
	for (size_t i = 0; i < A.rows_num; i++)
	{
		D.p[i][i] = A.p[i][i];
	}


	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix g0 = A * s0 - b;
	Matrix B0 = Matrix::eye(s0.rows_num);//初始时可以为任意对称正定矩阵
	Matrix B;
	Matrix d0 = B0 * g0;
	d0 = Matrix(d0.rows_num, d0.cols_num) - d0;
	Matrix temps;
	Matrix tempg;
	Matrix p;
	Matrix q;
	for (int k = 0; k < iteration; k++)
	{
		/*if (g0.norm2(0)<=eps)
		{
			break;
		}
		else*/
		{
			double alpha = -(Matrix::T(g0) * d0).p[0][0] / ((Matrix::T(d0) * A * d0)).p[0][0];

			temps = s0;//save s0

			s0 = s0 + d0 * alpha;//迭代

			tempg = g0;
			g0 = A * s0 - b;

			//更新B
			p = s0 - temps;
			q = g0 - tempg;
			//double gamma = (Matrix::T(q) * p).p[0][0] / (Matrix::T(q) * q).p[0][0];
			B = Matrix::inv(D);
			//使用双端队列来保存pi和qi
			deque_of_p.push_back(p);
			deque_of_q.push_back(q);
			deque_of_ro.push_back(1.0 / (Matrix::T(p) * q).p[0][0]);
			if (k > m)
			{
				//k>m删除前端元素
				deque_of_p.pop_front();
				deque_of_q.pop_front();
				deque_of_ro.pop_front();
			}
			d0 = LBFGS_two_loop(m, k, deque_of_p, deque_of_q
				, deque_of_ro, g0, B);
		}
	}
	return s0;
}
complex_matrix SD(complex_matrix H_hat,complex_matrix y_hat,double s,int iteration)//传统最速梯度下降法 使用梯度下降法求最小值
{
	complex_matrix y_1 = (!H_hat) * y_hat;
	complex_matrix A = ((!H_hat) * H_hat + complex_matrix::eyes(H_hat.cols_num) * s);
	double eps = 1e-5;
	complex_matrix  s0 = complex_matrix(H_hat.cols_num, 1);//generate_signal(H_hat.cols_num);//产生随机的s0
	for (size_t i = 0; i < iteration; i++)
	{
		s0 = s0 + y_1*(((!y_1)*y_1).p[0][0] / ((!(A*y_1))*y_1).p[0][0]);
		y_1 = y_1 - A * y_1 * (((!y_1) * y_1).p[0][0] / (((!(A * y_1))) * y_1).p[0][0]);
	}
	//s0.Show();
	return s0;
}


complex_matrix BB(complex_matrix H_hat, complex_matrix y_hat, double s,int iteration)// bb算法
{
	double eps = 1e-5;
	complex_matrix y_1 = (!H_hat) * y_hat;
	complex_matrix s0(H_hat.cols_num, 1);
	complex_matrix s1 =complex_matrix(s0.rows_num,1);
	for (int i = 0; i < s0.rows_num;i++)
	{
		s1.p[i][0] = Complex(1, 0);
	}
	complex_matrix x0 = s1 - s0;
	complex_matrix A = (!H_hat) * H_hat + complex_matrix::eyes(H_hat.cols_num) * s;

	complex_matrix temps=s1;
	for (size_t i = 0; i < iteration; i++)
	{
		if (x0.norm2(0) <= eps&&i!=0)
		{
			break;
		}
		else
		{
			temps = s1;
			s1 = s1 -gradient_complex(A,s1,y_1)*(x0.cdot() / (A*x0).cdot(x0));//已经为负梯度
			x0 = s1 - temps;
		}
	}
	return s1;
}
complex_matrix CG(complex_matrix H_hat, complex_matrix y_hat, double s,int iteration)//共轭梯度算法
{
	double eps = 1e-5;
	complex_matrix y_1 = (!H_hat) * y_hat;
	complex_matrix A = ((!H_hat) * H_hat + complex_matrix::eyes(H_hat.cols_num) * s);
	complex_matrix s0(H_hat.cols_num, 1);
	complex_matrix y0 = y_1;
	complex_matrix p0 = y0;
	complex_matrix tempy = y0;;
	for (size_t i = 0; i < iteration; i++)
	{
			tempy = y0;
			s0 = s0 + p0 * (((!y0) * y0).p[0][0] / ((!(A * p0)) * p0).p[0][0]);
			y0 = y0 - A * p0 * (((!y0) * y0).p[0][0] / ((!(A * p0)) * p0).p[0][0]);
			p0 = y0 + p0 * ((((!y0) * y0).p[0][0]) / ((!tempy) * tempy).p[0][0]);
	}
	return s0;
}
complex_matrix SPCG(complex_matrix H_hat, complex_matrix y_hat, double s,int iteration)//SPCG算法
{
	double eps = 1e-5;
	complex_matrix y_1 = (!H_hat) * y_hat;
	complex_matrix A = ((!H_hat) * H_hat + complex_matrix::eyes(H_hat.cols_num) * s);
	complex_matrix s0(H_hat.cols_num, 1);
	complex_matrix B(A.rows_num, A.cols_num);
	for (size_t i = 0; i < B.rows_num; i++)
	{
		B.p[i][i] =Complex( 1.0 / sqrt(A.p[i][i].real),0);//HT*H虚部为零0
	}

	complex_matrix y0 = B*y_1;
	complex_matrix p0 = y0;
	complex_matrix tempy = y0;
	for (size_t i = 0; i < iteration; i++)
	{
		tempy = y0;
		s0 = s0 + B * p0 * (y0.cdot() / (B * A * B * p0).cdot(p0));
		y0 = y0 - B * A * B * p0 * (y0.cdot() / (B * A * B * p0).cdot(p0));
		p0 = y0 + p0 * (y0.cdot() / tempy.cdot()); 
	}
	return s0;
}
Matrix Jacobi(complex_matrix H_hat, complex_matrix y_hat, double s=1.0, int iteration=1)//雅可比迭代
{
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;


	Matrix D(A.rows_num, A.cols_num);
	Matrix N_1(A.rows_num, A.cols_num); 
	for (size_t i = 0; i <A.rows_num; i++)
	{
		D.p[i][i] = A.p[i][i];
	}
	N_1 = D - A;
	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix temps = s0;
	for (size_t i = 0; i < iteration; i++)
	{
		/*if ((temps - s0).norm2(0) <= eps && i != 0)
		{
			break;
		}
		else*/
		{
			temps = s0;
			s0 = Matrix::inv(D) * (b + N_1 * s0);
		}
	}
	return s0;
}

Matrix Gauss(complex_matrix H_hat, complex_matrix y_hat, double s, int iteration = 1)//高斯赛德尔迭代
{
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;

	Matrix D(A.rows_num, A.cols_num);
	Matrix L(A.rows_num, A.cols_num);
	Matrix U(A.rows_num, A.cols_num);
	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix temps = s0;
	for (size_t i = 0; i < A.rows_num; i++)
	{
		D.p[i][i] = A.p[i][i];
	}
	for (size_t i = 0; i < A.rows_num; i++)
	{
		for (size_t j = 0; j < A.cols_num; j++)
		{
			if (i > j)
			{
				L.p[i][j] = A.p[i][j];
			}
			else if (i < j)
			{
				U.p[i][j] = A.p[i][j];
			}
		}
	}
	for (size_t i = 0; i < iteration; i++)
	{
		/*if ((temps - s0).norm2(0) <= eps && i != 0)
		{
			break;
		}
		else*/
		{
			temps = s0;
			s0 = Matrix::inv(D + L) * (b - U * s0);
		}
	}
	return s0;

}
Matrix SOR(complex_matrix H_hat, complex_matrix y_hat, double s, double omega, int iteration = 1)//超松弛迭代
{
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;

	Matrix D(A.rows_num, A.cols_num);
	Matrix L(A.rows_num, A.cols_num);
	Matrix U(A.rows_num, A.cols_num);
	Matrix s0 = generate_signal(K).realmatrix();//产生一个随机的s0
	Matrix temps = s0;
	for (size_t i = 0; i < A.rows_num; i++)
	{
		D.p[i][i] = A.p[i][i];
	}
	for (size_t i = 0; i < A.rows_num; i++)
	{
		for (size_t j = 0; j < A.cols_num; j++)
		{
			if (i > j)
			{
				L.p[i][j] = A.p[i][j];
			}
			else if (i < j)
			{
				U.p[i][j] = A.p[i][j];
			}
		}
	}
	for (size_t i = 0; i < iteration; i++)
	{
		/*if ((temps - s0).norm2(0) <= eps && i != 0)
		{
			break;
		}
		else*/
		{
			temps = s0;
			s0 = Matrix::inv(D + L * omega) * (D * (1 - omega) - U * omega) * s0 + Matrix::inv(D + L * omega) * b * omega;
		}
	}
	return s0;
}
Matrix Richard(complex_matrix H_hat, complex_matrix y_hat, double s, double omega, int iteration = 1)
{
	double eps = 1e-6;
	int N = 2 * y_hat.rows_num;
	int K = 2 * H_hat.cols_num;
	//complex_matrix temp_y(2 * y.rows_num, 1);
	Matrix R_y = y_hat.realmatrix();//实部矩阵
	Matrix I_y = y_hat.imagmatrix();//虚部矩阵

	Matrix R_H = H_hat.realmatrix();
	Matrix I_H = H_hat.imagmatrix();

	Matrix y = R_y.combine_rows(R_y, I_y);

	Matrix H = R_H.combine_rows(R_H, I_H);
	Matrix temp_imagpart = Matrix(H_hat.rows_num, H_hat.cols_num) - I_H;//负的H的虚部矩阵
	Matrix temp = temp_imagpart.combine_rows(temp_imagpart, R_H);
	H = H.combine_cols(H, temp);//最终得到H
	Matrix A = Matrix::T(H) * H + Matrix::eye(K) * s;
	Matrix b = Matrix::T(H) * y;

	Matrix s0 (K,1);//0矩阵
	Matrix temps = s0;
	
	for (size_t i = 0; i < 1000; i++)
	{
		if ((temps - s0).norm2(0) <= eps && i != 0)
		{
			break;
		}
		else
		{
			temps = s0;
			s0 = s0 + (b - A * s0) * omega;
		}
	}
	return s0;
}



