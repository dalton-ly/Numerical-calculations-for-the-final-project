#pragma once
#include"complex_matrix.h"
#include"auxiliary.h"
#include<fstream>
//#include<thread>
#include<iomanip>
/*****************直接求伪逆(ZF)*************************/
complex_matrix ZF(const complex_matrix& x, const complex_matrix& H)
{
	complex_matrix temp_x = x;
	complex_matrix c(H.cols_num, 1);
	complex_matrix inv_of_H = H.pseudo_inverse();//pseudo_inverse of H
	c = inv_of_H * x;
	for (int i = 0; i < c.rows_num; i++)
		c.p[i][0] = quantization(c.p[i][0]);
	return c;
}
//V_BLAST算法
complex_matrix V_BLAST(complex_matrix& x, complex_matrix& H)
{
	complex_matrix temp_x = x;
	complex_matrix temp_H = H;
	complex_matrix re(H.cols_num, 1);
	//int* flag = new int[H.rows_num];
	Complex c(0, 0);

	int** preserve_position = new int* [2];
	preserve_position[0] = new int[H.cols_num];
	preserve_position[1] = new int[H.cols_num];
	for (int i = 0; i < H.cols_num; i++)
	{
		preserve_position[0][i] = i;	//第一行记录列下标
		preserve_position[1][i] = 0;	//第二行记录当前列是否除去，初始为0(未处理)
	}
	int k = 0;
	for (int i = 0; i < H.cols_num; i++)
	{
		int* original_position = new int[H.cols_num - i];
		k = 0;
		for (int j = 0; j < H.cols_num; j++)
		{
			if (!preserve_position[1][j])//如果这列没有被除去
			{
				original_position[k] = preserve_position[0][j];//记录在原矩阵的位置 
				k++;
			}
		}
		complex_matrix G = temp_H.pseudo_inverse();
		//找出G中2-范数最小的一行
		int min_row = 0;//从第一行开始
		complex_matrix g = G.get_row(0);
		double norm_min = g.norm2_vector();//第一行的2范数
		double temp_norm;
		for (int j = 1; j < G.rows_num; j++)
		{
			g = G.get_row(j);
			temp_norm = g.norm2_vector();
			if (temp_norm < norm_min)
			{
				min_row = j;
				norm_min = temp_norm;
			}
		}
		complex_matrix w_T = G.get_row(min_row);
		complex_matrix y = w_T * temp_x;
		c = quantization(y.p[0][0]);
		//original_position [min_row]的值为min_row在原H矩阵的位置
		re.p[original_position[min_row]][0] = c;
		temp_x = temp_x - (temp_H.get_col(min_row)) * c;
		temp_H = temp_H.remove_col(min_row);
		//标记min-row行已经被处理
		preserve_position[1][original_position[min_row]] = 1;
		
		//preserve_position[0][original_position[min_row]] = i;
		//释放内存
		delete[]original_position;
	}
	for (int i = 0; i < temp_x.rows_num; i++)
		temp_x.p[i][0] = quantization(temp_x.p[i][0]);
	return re;
}
complex_matrix QRD(complex_matrix& x, complex_matrix& H)
{
	complex_matrix c_hat(H.cols_num, 1);		//还原得到的信号
	complex_matrix H_temp = H;	//H的复制，防止修改原本H
	complex_matrix y;		//接收信号x被Q修正后的信号
	complex_matrix Q_H;		//Q的共轭转置
	//QR分解
	complex_matrix Q = complex_matrix::gram_improve(H_temp)[0];
	complex_matrix R = complex_matrix::gram_improve(H_temp)[1];
	//x根据Q变换为y
	Q_H = !Q;
	y = Q_H * x;
	//利用R的nT行和y还原出发射信号
	//向前回代还原
	for (int k = c_hat.rows_num - 1; k >= 0; k--)
	{
		Complex dk;
		Complex zk;
		for (int i = k + 1; i < H.cols_num; i++)
		{
			dk = dk + R.p[k][i] * c_hat.p[i][0];
		}
		zk = y.p[k][0] - dk;
		c_hat.p[k][0] = quantization(zk / R.p[k][k]);
	}
	return c_hat;
}
/**********选择QR*******/
complex_matrix SQRD(complex_matrix& x, complex_matrix& H)
{
	int* S = new int[H.cols_num];//S记录交换顺序
	for (int r = 0; r < H.cols_num; r++ )
	{
		S[r] = r;
	}
	//保存数据防止修改
	complex_matrix temp_x = x;
	complex_matrix temp_H = H;
	complex_matrix c_hat(H.cols_num, 1);
	//SQRD算法部分
	complex_matrix R(H.cols_num,H.cols_num);
	complex_matrix Q = H;
	for (int i = 0; i < H.cols_num; i++)
	{
		//查找范数最小列的下标
		complex_matrix ql = Q.get_col(i);
		double min_norm = ql.norm2_vector();
		int k_i = i;
		for (int l = i+1; l < H.cols_num; l++)
		{
			complex_matrix temp_col_vector = Q.get_col(l);
			double temp_norm = temp_col_vector.norm2_vector();
			if (temp_norm < min_norm)
			{
				k_i = l;
				min_norm = temp_norm;
			}
		}
		//exchange Q R S
		R.exchange_col(i, k_i);
		Q.exchange_col(i, k_i);
		int S_change_temp = S[i];
		S[i] = S[k_i];
		S[k_i] = S_change_temp;

		double r_ii = Q.get_col(i).norm2_vector();
		R.p[i][i].real = r_ii;
		//qi=qi/rii
		for (int t = 0; t < Q.rows_num; t++)
		{
			Q.p[t][i]=Q.p[t][i] / r_ii;
		}
		//
		Complex r_il;
		for (int l = i + 1; l < H.cols_num; l++)
		{
			r_il =( !(Q.get_col(i)) * Q.get_col(l)).p[0][0];
			R.p[i][l] = r_il;
			//更改ql
			for (int t = 0; t < Q.rows_num; t++)
			{
				Complex temp=Q.p[t][i] * r_il;
				Q.p[t][l] = Q.p[t][l] - temp;
			}
		}
	}
	//信号检测部分
	complex_matrix y = (!Q) * temp_x;
	for (int k = H.cols_num - 1; k >= 0; k--)
	{
		//求dk
		Complex dk;
		Complex zk;
		for (int i = k + 1; i < H.cols_num; i++)
		{
			dk = dk + R.p[k][i] * c_hat.p[i][0];
		}
		zk = y.p[k][0] - dk;
		c_hat.p[k][0] = quantization(zk / R.p[k][k]);
	}
	//根据S还原C_hat
	complex_matrix c_temp = c_hat;
	for (int a = 0; a < c_hat.rows_num; a++)
	{
		c_hat.p[S[a]][0] = c_temp.p[a][0];
		//c_hat.p[a][0] = c_temp.p[S[a]][0];  
	}
	return c_hat;
}
complex_matrix MMSE(complex_matrix& x, complex_matrix& H,double S)//S为噪声标准差
{
	complex_matrix temp_x = x;
	complex_matrix temp_H = H;
	complex_matrix s_hat(H.cols_num, 1);//结果矩阵
	complex_matrix x0(H.cols_num, 1);//nT*1的0向量
	complex_matrix I = (complex_matrix::eyes(H.cols_num)*S);
	complex_matrix H_ = temp_H.combine_rows(temp_H, I);
	complex_matrix x_ = temp_x.combine_rows(temp_x, x0);
	s_hat = H_.pseudo_inverse() * x_;
	for (int i = 0; i < H.cols_num; i++)//进行量化
	{
		s_hat.p[i][0] = quantization(s_hat.p[i][0]);
	}
	return s_hat;
}
//无排序MMSEQR
complex_matrix MMSE_QR(complex_matrix & x, complex_matrix & H,double S)
{
	complex_matrix temp_x = x;
	complex_matrix temp_H = H;
	complex_matrix x0(H.cols_num, 1);
	complex_matrix c_hat(H.cols_num, 1);//结果矩阵
	complex_matrix I = complex_matrix::eyes(H.cols_num)*S;
	complex_matrix H_ = temp_H.combine_rows(temp_H, I);
	complex_matrix x_ = temp_x.combine_rows(temp_x, x0);
	complex_matrix Q_ = complex_matrix::gram_improve(H_)[0];
	complex_matrix R_= complex_matrix::gram_improve(H_)[1];
	complex_matrix y = (!Q_) * x_;
	for (int k = c_hat.rows_num - 1; k >= 0; k--)
	{
		Complex dk;
		Complex zk;
		for (int i = k + 1; i < H.cols_num; i++)
		{
			dk = dk + R_.p[k][i] * c_hat.p[i][0];
		}
		zk = y.p[k][0] - dk;
		c_hat.p[k][0] = quantization(zk / R_.p[k][k]);
	}
	return c_hat;
}
complex_matrix ML(complex_matrix& x, complex_matrix& H)//ML检测算法
{
	complex_matrix temp_x = x;
	complex_matrix temp(H.cols_num, 1);
	int* bits = new int[2 * H.cols_num];
	double min_norm = temp_x.norm2_vector();
	int k = 0;//记录对应的01序列的10进制值
	//遍历星座图上所有可能的发射向量，查找最小欧式距离对应的发射向量
	for (int i = pow(2, H.cols_num * 2)-1; i >=0; i--)
	{
		//计算||x-Hx_1|| x_1为所有可能的情况
		int a = i;
		for (int j = 0; j < H.cols_num*2; j++)//zhuan hua wei er jin zhi
		{
			bits[j] = a % 2;
			a = a / 2;
		}
		//将bits转化为复数矩阵
		for (int t = 0; t < H.cols_num * 2; t++)
		{
			temp.p[t / 2][0].real = bits[t] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
			t++;
			//t+1/2由于整形仍为t/2
			temp.p[t / 2][0].imag = bits[t] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
		}
		double m = (temp_x-H*temp).norm2_vector();
		if (m < min_norm)
		{
			min_norm = m;
			k = i;
		}
	}
	for (int i = 0; i < H.cols_num*2; i++)//2norm最小的x
	{
		bits[i] = k% 2;
		k = k / 2;
	}
	for (int t = 0; t < H.cols_num * 2; t++)
	{
		temp.p[t / 2][0].real = bits[t] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
		t++;
		//t+1/2由于整形仍为t/2
		temp.p[t / 2][0].imag = bits[t] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
	}
	return temp;
}

//ofstream output;
//void writeBERtofile(int nR, int nT, int times)
//{
//	//time 为测试次数
//	//对于7种算法分别测试times次每次改变SNR信噪比，每次测试得到的BER保存在一个数组中，测试结束后求数组的平均值
//	//然后写入到文件中，文件名为方法名称对应的BER
//	//output.open("test.txt");
//	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
//	int SNR_max = 13;
//	//以下向量用于保存每个算法的平均误码率，随后写入文件
//	Matrix average_of_VBLAST_vs_SNR(SNR_max, 1);
//	Matrix average_of_QRD_vs_SNR(SNR_max, 1);
//	Matrix average_of_SQRD_vs_SNR(SNR_max, 1);
//	Matrix average_of_ZF_vs_SNR(SNR_max, 1);
//	Matrix average_of_MMSE_vs_SNR(SNR_max, 1);
//	Matrix average_of_MMSE_QR_vs_SNR(SNR_max, 1);
//	Matrix average_of_ML_vs_SNR(SNR_max, 1);
//	double* SNR_1 = new double[SNR_max];
//
//	for (int i = 0; i < SNR_max; i++)
//		SNR_1[i] = SNR + i;
//	complex_matrix* all_H = new complex_matrix[times];//所有用于计算的H
//	int** all_bits = new int* [times];//所有比特信号
//	srand((unsigned)time(NULL));
//	for (int t = 0; t < times; t++)
//	{
//		all_H[t] = generate_H(nR,nT);
//		all_bits[t] = generate_random_bits(nT);
//	}
//
////#pragma omp parallel for
//		for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
//		{
//			//double Eb_N0 = pow(10, SNR / 10);
//			double BER_of_VBLAST = 0;
//			double BER_of_QRD = 0;
//			double BER_of_SQRD = 0;
//			double BER_of_MMSE = 0;
//			double BER_of_ZF = 0;
//			double BER_of_MMSE_QR = 0;
//			double BER_of_ML = 0;
//			srand((unsigned)time(NULL));
//			for (int i = 0; i < times; i++)
//			{
//				//考虑使用多线程，同时计算七个结果
//				complex_matrix temp_H    = all_H[i] * sqrt(SNR_1[k]);//设置H矩阵的方差
//				complex_matrix c         = generate_signal(nT);//发射信号
//				complex_matrix miu       = generate_noise(nR);
//				complex_matrix temp_x = temp_H * c + miu;
//				complex_matrix result_of_VBLAST = V_BLAST(temp_x, temp_H);
//				complex_matrix result_of_QRD = QRD(temp_x, temp_H);
//				complex_matrix result_of_SQRD = SQRD(temp_x, temp_H);
//				complex_matrix result_of_ZF = ZF(temp_x, temp_H);
//				complex_matrix result_of_MMSE = MMSE(temp_x, temp_H);
//				complex_matrix result_of_MMSE_QR = MMSE_QR(temp_x, temp_H);
//				complex_matrix result_of_ML = ML(temp_x, temp_H);
//
//				BER_of_VBLAST += BitsErrorRate(c, result_of_VBLAST);
//				BER_of_QRD += BitsErrorRate(c, result_of_QRD);
//				BER_of_SQRD += BitsErrorRate(c, result_of_SQRD);
//				BER_of_MMSE += BitsErrorRate(c, result_of_MMSE);
//				BER_of_ZF += BitsErrorRate(c, result_of_ZF);
//				BER_of_MMSE_QR += BitsErrorRate(c, result_of_MMSE_QR);
//				BER_of_ML += BitsErrorRate(c, result_of_ML);
//			}
//			average_of_VBLAST_vs_SNR.setpoint(k, 0, BER_of_VBLAST / times);
//			average_of_QRD_vs_SNR.setpoint(k, 0, BER_of_QRD / times);
//			average_of_SQRD_vs_SNR.setpoint(k, 0, BER_of_SQRD / times);
//			average_of_ZF_vs_SNR.setpoint(k, 0, BER_of_ZF / times);
//			average_of_MMSE_vs_SNR.setpoint(k, 0, BER_of_MMSE / times);
//			average_of_MMSE_QR_vs_SNR.setpoint(k, 0, BER_of_MMSE_QR / times);
//			average_of_ML_vs_SNR.setpoint(k, 0, BER_of_ML / times);
//		}
//	
//	output.open("BER_of_VBLASR.txt");//,ios::app);
//	for (int a = 0; a < SNR_max; a++)
//	{
//		output << setprecision(10) << average_of_VBLAST_vs_SNR.Point(a, 0) << "  ";
//	}
//	output.close();
//	output.open("BER_of_QRD.txt");//,ios::app);
//	for (int a = 0; a < SNR_max; a++)
//	{
//		output << setprecision(10) << average_of_QRD_vs_SNR.Point(a, 0) << "  ";
//	}
//	output.close();
//	output.open("BER_of_SQRD.txt");//,ios::app);
//	for (int a = 0; a < SNR_max; a++)
//	{
//		output << setprecision(10) << average_of_SQRD_vs_SNR.Point(a, 0) << "  ";
//	}output.close();
//	output.open("BER_of_MMSE.txt");//,ios::app);
//	for (int a = 0; a < SNR_max; a++)
//	{
//		output << setprecision(10) << average_of_MMSE_vs_SNR.Point(a, 0) << "  ";
//	}output.close();
//	output.open("BER_of_MMSE_QR.txt");//,ios::app);
//	for (int a = 0; a < SNR_max; a++)
//	{
//		output <<setprecision(10)<< average_of_MMSE_QR_vs_SNR.Point(a, 0) << "  ";
//	}output.close();
//	output.open("BER_of_ZF.txt");//,ios::app);
//	for (int a = 0; a < SNR_max; a++)
//	{
//		output << setprecision(10) << average_of_ZF_vs_SNR.Point(a, 0) << "  ";
//	}output.close();
//output.open("BER_of_ML.txt");//,ios::app);
//for (int a = 0; a < SNR_max; a++)
//	{
//		output << setprecision(10) << average_of_ML_vs_SNR.Point(a, 0) << "  ";
//	}output.close();
//}
//使用多线程，同时计算七个结果
ofstream outVBLAST;
void VBLAST_thread(int nR, int nT, int times)
{
	outVBLAST.open("1_BER_of_VBLAST_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outVBLAST << "Eb_N0" << "			" << "BER" << endl;
	for (int k =0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_VBLAST = 0;;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu *(2)/sqrt(pow(10,double(k)/10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_VBLAST = V_BLAST(temp_x, temp_H);
			BER_of_VBLAST += BitsErrorRate(c, result_of_VBLAST);
		} 
		outVBLAST << fixed << setprecision(10) << k<<"			" << BER_of_VBLAST / times << "  " << endl;;
	}
	outVBLAST.close();
}
ofstream outSQRD;
void SQRD_thread(int nR, int nT, int times)
{
	outSQRD.open("1_BER_of_SQRD_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outSQRD << "Eb_N0" << "			" << "BER" << endl;
	for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_SQRD = 0;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu*2 / sqrt(pow(10, double(k) /10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_SQRD = SQRD(temp_x, temp_H);
			BER_of_SQRD += BitsErrorRate(c, result_of_SQRD);
		}
		outSQRD << fixed << setprecision(10) << k <<"			"<< BER_of_SQRD / times << "  " << endl;;
	}
	outSQRD.close();
}
ofstream outQRD;
void QRD_thread(int nR, int nT, int times)
{
	outQRD.open("1_BER_of_QRD_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outQRD << "Eb_N0" << "			" << "BER" << endl;
	for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_QRD = 0.0;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu*2 / sqrt(pow(10, double(k) / 10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_QRD = QRD(temp_x, temp_H);
			BER_of_QRD += BitsErrorRate(c, result_of_QRD);
		}
		outQRD << fixed << setprecision(10) << k <<"			"<< BER_of_QRD / times << "  " << endl;;
	}
	outQRD.close();
}
ofstream outZF;
void ZF_thread(int nR, int nT, int times)
{
	outZF.open("1_BER_of_ZF_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outZF << "Eb_N0" << "			" << "BER" << endl;
	for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_ZF = 0;;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu*2 / sqrt(pow(10, double(k) / 10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_ZF = ZF(temp_x, temp_H);
			BER_of_ZF += BitsErrorRate(c, result_of_ZF);
		}
		outZF << fixed << setprecision(10) << k <<"			"<< BER_of_ZF / times << "  " << endl;;
	}
	outZF.close();
}
ofstream outMMSE;
void MMSE_thread(int nR, int nT, int times)
{
	outMMSE.open("1_BER_of_MMSE_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outMMSE << "Eb_N0" << "			" << "BER" << endl;
	for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_MMSE = 0;;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu*2 / sqrt(pow(10, double(k) / 10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_MMSE = MMSE(temp_x, temp_H, double(2) /sqrt(pow(10, double(k) / 10.0)));
			BER_of_MMSE += BitsErrorRate(c, result_of_MMSE);
		}
		outMMSE << fixed << setprecision(10) << k<<"			" << BER_of_MMSE / times << "  " << endl;;
	}
	outMMSE.close();
}
ofstream outMMSE_QR;
void MMSE_QR_thread(int nR, int nT, int times)
{
	outMMSE_QR.open("1_BER_of_MMSE_QR_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outMMSE_QR << "Eb_N0" << "	" << "BER" << endl;
	for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_MMSE_QR = 0;;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu *2/ sqrt(pow(10, double(k) / 10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_MMSE_QR = MMSE_QR(temp_x, temp_H, double(2)/ sqrt(pow(10, double(k)/10.0)));
			BER_of_MMSE_QR += BitsErrorRate(c, result_of_MMSE_QR);
		}
		outMMSE_QR << fixed << setprecision(10) << k<<"			" << BER_of_MMSE_QR / times << endl;;
	}
	outMMSE_QR.close();
}
ofstream outML;
void ML_thread(int nR, int nT, int times)
{
	outML.open("1_BER_of_ML_thread.txt");
	double SNR = 1.0;//信噪比 SNR=10lg(Eb/N0)
	int SNR_max = 13;
	srand((unsigned)time(NULL));
	complex_matrix temp_H = generate_H(nR, nT);
	outML << "Eb_N0" << "	" << "BER" << endl;
	for (int k = 0; k < SNR_max; k++)//Eb/N0等于20 SNR约为13
	{
		double BER_of_ML = 0;;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//发射信号
			for (int v = 0; v < c.rows_num; v++)
			{
				c.p[v][0] = quantization(c.p[v][0]);
			}
			complex_matrix miu = generate_noise(nR);
			miu = miu *2/ sqrt(pow(10, double(k) / 10.0));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result_of_ML = ML(temp_x, temp_H);
			BER_of_ML += BitsErrorRate(c, result_of_ML);
		}
		outML << fixed << setprecision(10) << k <<"			"<< BER_of_ML / times << "  " << endl;;
	}
	outML.close();
}













//总体误码率测试：指定发射和接收天线数（输入参数先接收，后发射），产生n个H信道
//每次用这n个H跑一个采样点，随后改变信噪比，再跑下一次same_SNR_n_H
//存储误码率采用一个行向量，分别记录每种方法，由一个信噪比和一种方法可确定该方法在此信噪比下的平均误码率
//采样点的数据都输出到txt文件中以便MATLAB画图处理
//一个采样点：指定发射和接收天线数，输入一个信噪比（功率比）
//n个H以该信噪比跑一次one_H_nn，返回n个H用每种方法时的平均误码率，即每种方法在该点的误码率BER
//一个采样点的一部分：指定发射和接收天线数（在H的维度信息里体现），同一个H
//每种方法“同步”跑nn次（源比特流和噪声都不同），返回每种方法的平均误码率

