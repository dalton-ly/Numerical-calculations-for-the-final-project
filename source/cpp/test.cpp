#include"algrithm_mimo.h"
#include"auxiliary.h"
#include<ctime>
#include<thread>
#include"ÆÚÄ©mimo.h"
#include<fstream>
int main()
{
	srand(time(NULL));
	int times = 1000;
	int nT = 16;
	int nR = 128;

	int SNR_max = 100;
	double BER_of_BFGS;
	complex_matrix temp_H = generate_H(nR, nT);
	for (int k = 1; k < SNR_max; k += 5)
	{
		srand((unsigned)time(NULL));
		double BER = 0;
		for (int i = 0; i < times; i++)
		{
			complex_matrix c = generate_signal(nT);//·¢ÉäÐÅºÅ
			c = quantization(c);
			complex_matrix miu = generate_noise(nR);
			miu = (miu / sqrt(k));
			complex_matrix temp_x = temp_H * c + miu;
			complex_matrix result = Matrix_to_complex(Jacobi(temp_H, temp_x, 0.5 * double(k), 5));
			BER += BitsErrorRate(c, result);
		}
		cout << BER << endl;

	}
}