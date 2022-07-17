#include"auxiliary.h"
int* generate_random_bits(int nT)//0 1 sequences
{
	double* random_nums = new double[2 * nT];
	int* result = new int[2 * nT];
	//srand((int)time(0));
	for (int i = 0; i < 2 * nT; i++)
	{
		random_nums[i] = ((double)rand() / RAND_MAX) * 2 - 1;//正负随机
		if (random_nums[i] >= 0)
			result[i] = 1;
		else
			result[i] = 0;
	}
	delete[]random_nums;
	return result;
}
complex_matrix generate_signal(int nT)//generate a complex matrix产生随机的矩阵
{
	complex_matrix signal(nT, 1);
	for (int i = 0; i < nT; i++)
	{
			signal.p[i][0].real = gaussianrand()/sqrt(2);
			signal.p[i][0].imag = gaussianrand()/sqrt(2);
	}
	return signal;
}
complex_matrix generate_signal(int* bits,int nT)//generate a complex matrix产生随机的矩阵
{
	complex_matrix signal(nT, 1);
	for (int i = 0; i < nT * 2; i++)
	{
		signal.p[i / 2][0].real = bits[i] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
		i++;
		signal.p[i / 2][0].imag = bits[i] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
	}
	return signal;
}

//产生高斯分布 box-muller法
double gaussianrand(double V)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	//srand((unsigned)time(NULL));
	if (phase == 0)
	{
		do {
			
			//获得两个（-1,1）的独立随机变量
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);	//由于涉及求对数，舍去S=0这个点情况，不影响连续随机变量的分布
		//Box - Muller转换
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;	//两个X都服从均值为0，方差为1的高斯分布，但每次选其中一个，下次选另一个
	//指定方差和期望值X = X * V + E;
	X = X * V ;
	return X;
}
//quantization 
Complex quantization(Complex y)
{
	Complex result;
	//第一象限和正半实轴上的估算为(1+j)/sqrt(2)
	if ((y.real > 0) && (y.imag >= 0)) 
	{
		result.real = 1.0 / sqrt(2.0); result.imag = 1.0 / sqrt(2.0); 
	}
	//第二象限和正半虚轴上的估算为(-1+j)/sqrt(2)
	else if ((y.real <= 0) && (y.imag > 0))
	{ 
		result.real = -1.0 / sqrt(2.0); result.imag = 1.0 / sqrt(2.0);
	}
	//第三象限和负半实轴上的估算为(-1-j)/sqrt(2)
	else if ((y.real < 0) && (y.imag <= 0))
	{ 
		result.real = -1.0 / sqrt(2.0); result.imag = -1.0 / sqrt(2.0);
	}
	//第四象限和负半虚轴上的估算为(1-j)/sqrt(2)
	else if ((y.real >= 0) && (y.imag < 0))
	{ 
		result.real = 1.0 / sqrt(2.0); result.imag = -1.0 / sqrt(2.0);
	}
	//原点上的估算为1+j
	else
	{
		result.real = 1.0 / sqrt(2.0); result.imag = 1.0 / sqrt(2.0);
	}
	return result;
}
int* signal_to_bits(complex_matrix& c)//将信号转化为比特
{
	int* bits = new int[c.rows_num * 2];
	for (int i = 0; i < c.rows_num * 2; i++)
	{
		//real<-1/sqrt(2) 1
		bits[i] = (abs(c.p[i / 2][0].real + 1 / sqrt(2.0)) < 1e-16) ? 1 : 0;
		i++;
		bits[i] = (abs(c.p[i / 2][0].imag + 1 / sqrt(2.0)) < 1e-16) ? 1 : 0;
	}
	return bits;
}
//产生一个高斯随机变量信道矩阵（均值0，方差默认1）
complex_matrix generate_H(int nR, int nT)
{
	complex_matrix H(nR, nT);
	for (int i = 0; i < H.rows_num; i++)
	{
		for (int j = 0; j < H.cols_num; j++)
		{
			H.p[i][j].real = gaussianrand()/sqrt(2);
			H.p[i][j].imag = gaussianrand()/sqrt(2);
		}
	}
	return H;
}

//产生一个高斯随机变量噪声向量（均值0，方差默认为1）
complex_matrix generate_noise(int nR)
{
	complex_matrix miu(nR, 1);
	for (int i = 0; i < miu.rows_num; i++)
	{
		miu.p[i][0].real = gaussianrand()/ sqrt(2);
		miu.p[i][0].imag = gaussianrand()/sqrt(2);
	}
	return miu;
}
/*complex_matrix rand_matrix(const complex_matrix& A)
{
	
	complex_matrix temp = A;
	for (int i = 0; i < A.rows_num; i++)
	{
		for (int j = 0; j < A.cols_num; j++)
		{
			temp.p[i][j].real = gaussianrand();
			temp.p[i][j].imag = gaussianrand();
		}
	}
	return temp;
}*/
//误码率分析
//比较两个比特流，得出误码率
double BitsErrorRate(complex_matrix origin_signal , complex_matrix signal_hat)
{
	////量化
	//for (int u = 0; u < origin_signal.rows_num;u++)
	//{
	//	origin_signal.p[u][0] = quantization(origin_signal.p[u][0]);
	//}
	signal_hat = quantization(signal_hat);
	int* bits = signal_to_bits(origin_signal);
	int* bits_hat = signal_to_bits(signal_hat);
	int nT = origin_signal.rows_num;
	int error = 0;
	for (int i = 0; i < nT * 2; i++)
	{
		if (bits[i] != bits_hat[i])
		{
			error++;
		}
	}
	double BER = double(error) / (double(nT) * 2.0);
	return BER;
}
double partial(double (*f)(Matrix, Matrix, Matrix ), Matrix x,Matrix A,Matrix b,int i)//求f函数在第i个变量的偏导
{
	double eps = 1e-8;
	Matrix temp=x;
	temp.p[i][0] = temp.p[i][0] + eps;
	double res = (*f)(temp, A,b) - (*f)(x, A,b);
	return res / eps;
}
Complex partial(Complex (*f)(complex_matrix, complex_matrix,
	complex_matrix,complex_matrix), complex_matrix x, complex_matrix A, complex_matrix b
	,complex_matrix y, int i)//求f函数在第i个变量的偏导
{
	Complex eps(1e-8,0);
	complex_matrix temp = x;
	temp.p[i][0] = temp.p[i][0] + eps;
	Complex res= (*f)(temp, A, b,y) - (*f)(x, A, b,y);
	return res/eps;
}
Matrix get_negetive_gradient(double (*f)(Matrix, Matrix,Matrix), Matrix x, Matrix A,Matrix b)//返回函数在一个x向量处的梯度矩阵 此作业中函数特指瑞丽商
{
	int rows = x.rows_num;
	Matrix gradient(rows, 1);

	for (int i = 0; i < rows; i++)
	{
		gradient.p[i][0] = -partial((*f), x, A, b,i);
	}
	return gradient;
}
complex_matrix get_negetive_gradient(Complex (*f)(complex_matrix, complex_matrix, complex_matrix,complex_matrix), complex_matrix x, complex_matrix A,
	complex_matrix b,complex_matrix y)//返回函数在一个x向量处的梯度矩阵 
{
	// 复数矩阵
	int rows = x.rows_num;
	complex_matrix gradient(rows, 1);

	for (int i = 0; i < rows; i++)
	{
		gradient.p[i][0] = Complex(0,0)-partial((*f), x, A, b,y, i);
	}
	return gradient;
}
complex_matrix gradient_complex(complex_matrix A, complex_matrix s, complex_matrix b)
{
	complex_matrix result;
	result = A * s - b;
	return result;
}
complex_matrix quantization(complex_matrix A)
{
	complex_matrix result(A.rows_num, A.cols_num);
	for (size_t i = 0; i < A.rows_num; i++)
	{	for (size_t j = 0; j < A.cols_num; j++)
		{
		result.p[i][j] = quantization(A.p[i][j]);
		}
	}
	return result;
}
complex_matrix Matrix_to_complex(Matrix A)
{
	complex_matrix temp(A.rows_num / 2,1);
	for (size_t i = 0; i < (A.rows_num/2); i++)
	{
		temp.p[i][0].real = A.p[i][0];
		temp.p[i][0].imag = A.p[i+A.rows_num/2][0];
	}
	return temp;
}