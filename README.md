# 11.25.2
LU
// homework1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<math.h>
#include<string>
#include <iomanip>
#include<fstream>
using namespace std;
double C[5][501];//压缩存储实对称矩阵A
const double  b = 0.16, c = -0.064;
const double epsilon = 1.0e-12;
void reverse(const double b, const double c);//求压缩矩阵
double PowerMethod(double lambda);//幂法.lambda为位移量
double InvPowerMethod(double lambda);//反幂法
double* Doolittle(double *B, double lambda);//Doolittle分解
double DoolittleDet();//LU分解求矩阵A的行列式的值
int max(int x, int y, int z);//求三个整数最大值
int min(int x, int y);
int max(int x, int y);
int main()
{
	reverse( b, c);
	double lambdaMax1, lambdaMax2;
	double lambda1, lambda501, lambdaS,lambdaik;
	lambdaMax1 = PowerMethod(0);
	lambdaMax2 = PowerMethod(lambdaMax1);
	if (lambdaMax1 > 0) {
		lambda501 = lambdaMax1 ;
		lambda1 = lambdaMax1+ lambdaMax2;
	}
	else {
		lambda1 = lambdaMax1;
		lambda501 = lambdaMax1 + lambdaMax2;
	}
	cout << setiosflags(ios::scientific)<< setprecision(12)<< "λ1=" << lambda1 << "， λ501=" << lambda501 << endl;
	lambdaS = InvPowerMethod(0);
	cout << "λs=" << lambdaS << endl;
	cout << "det(A)=" << DoolittleDet() << endl;
	cout << "cond(A)2=" << lambdaMax1 / lambdaS << endl;
	int changeline = 0;
	for (int k = 1; k <= 39; k++) {
		lambdaik = InvPowerMethod(lambda1+k*(lambda501 - lambda1) / 40);
		cout << "λi" << k << "=" << lambdaik<<"  ";
		changeline++;
		if (changeline % 4 == 0)
			cout << endl;
	}
	//ofstream write; //write只是个名字 你可以定义为任何其他的名字
	//write.open("text.txt"); //表示你要把内容输出到“text.txt"这个文件里 如果没有这个文件，会自动创建这个文件
	//write << setiosflags(ios::scientific) << setprecision(12) << "λ1=" << lambda1 << "， λ501=" << lambda501 << endl;
	//lambdaS = InvPowerMethod(0);
	//write << "λs=" << lambdaS << endl;
	//write << "det(A)=" << DoolittleDet() << endl;
	//write << "cond(A)2=" << lambdaMax1 / lambdaS << endl;
	//int changeline = 0;
	//for (int k = 1; k <= 39; k++) {
	//	lambdaik = InvPowerMethod(lambda1 + k*(lambda501 - lambda1) / 40);
	//	write << "λi" << k << "=" << lambdaik << "  ";
	//	changeline++;
	//	if (changeline % 4 == 0)
	//		write << endl;
	//}
    return 0;
}

void reverse(const double b, const double c) {
	double a[501];
	for (int i = 0; i < 501; i++) {
		a[i] = (1.64 - 0.024*(i+1))*sin(0.2*(i + 1)) - 0.64*exp(0.1 / (i + 1));
	}
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 501; j++) {
			switch (i)
			{
			case(0):
				if (j < 2) C[i][j] = 0;
				else C[i][j] = c;
				break;
			case(1):
				if (j == 0) C[i][j] = 0;
				else C[i][j] = b;
				break;
			case(2):
				C[i][j] = a[j];
				break;
			case(3):
				if (j < 500) C[i][j] = b;
				else C[i][j] = 0;
				break;
			case(4):
				if (j < 499) C[i][j] = c;
				else C[i][j] = 0;
				break;
			default:
				break;
			}
		}
	}
}
double PowerMethod(double lambda) {
	double u1[501],u2[501];
	for (int i = 0; i < 501; i++) {
		u1[i] = 1;
		u2[i] = 0;
	}
	u1[0] = 1.0;
	double eta=0, y[501], beta1=0,beta2=0;
	for (int i = 0; i < 501; i++) {//求二范数
		eta += u1[i]*u1[i];
	}
	eta = sqrt(eta);
	double CC[5][501];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 501; j++) {
			if (i == 2) {
				CC[i][j] = C[i][j] - lambda;
			}
			else {
				CC[i][j] = C[i][j];
			}
		}
	}
	for (int i = 0; i < 501; i++) {//归一化
		y[i] = u1[i] / eta;
	}
	for (int i = 0; i < 501; i++) {
		u2[i] = CC[2][i] * y[i] + CC[1][i] * y[i - 1] + CC[3][i] * y[i + 1] + CC[0][i] * y[i - 2] + CC[4][i] * y[i + 2];
	}
	for (int i = 0; i < 501; i++) {
		beta1 += u2[i] * y[i];
	}
	beta2 = beta1;
	do {
		beta1 = beta2;
		beta2 = 0;//beta2置零
		for (int i = 0; i < 501; i++) {//uk变为uk-1
			u1[i] = u2[i];
		}
		for (int i = 0; i < 501; i++) {//uk置零
			u2[i] = 0;
		}
		eta = 0;//eta置零
		for (int i = 0; i < 501; i++) {//求二范数
			eta += u1[i] * u1[i];
		}
		eta = sqrt(eta);
		for (int i = 0; i < 501; i++) {//归一化
			y[i] = u1[i] / eta;
		}
		for (int i = 0; i < 501; i++) {
			u2[i] = CC[2][i] * y[i] + CC[1][i] * y[i - 1] + CC[3][i] * y[i + 1] + CC[0][i] * y[i - 2] + CC[4][i] * y[i + 2];
		}
		for (int i = 0; i < 501; i++) {
			beta2 += u2[i] * y[i];
		}
	} while ((abs(beta2 - beta1) / abs(beta2)) > epsilon);
	return beta2;
}

double InvPowerMethod(double lambda) {
	double u1[501], u2[501];
	for (int i = 0; i < 501; i++) {
		u1[i] = 1;
		u2[i] = 0;
	}
	u1[0] = 1.0;
	double eta = 0, y[501], beta1 = 0, beta2 = 0;
	for (int i = 0; i < 501; i++) {//求二范数
		eta += u1[i] * u1[i];
	}
	eta = sqrt(eta);
	for (int i = 0; i < 501; i++) {//归一化
		y[i] = u1[i] / eta;
	}
	for (int i = 0; i < 501; i++) {
		u2[i] = Doolittle(y,lambda)[i];
	}
	for (int i = 0; i < 501; i++) {
		beta1 += u2[i] * y[i];
	}
	beta2 = beta1;
	do {
		beta1 = beta2;
		beta2 = 0;//beta2置零
		for (int i = 0; i < 501; i++) {//uk变为uk-1
			u1[i] = u2[i];
		}
		for (int i = 0; i < 501; i++) {//uk置零
			u2[i] = 0;
		}
		eta = 0;//eta置零
		for (int i = 0; i < 501; i++) {//求二范数
			eta += u1[i] * u1[i];
		}
		eta = sqrt(eta);
		for (int i = 0; i < 501; i++) {//归一化
			y[i] = u1[i] / eta;
		}
		for (int i = 0; i < 501; i++) {
			u2[i] = Doolittle(y, lambda)[i];
		}
		for (int i = 0; i < 501; i++) {
			beta2 += u2[i] * y[i];
		}
	} while ((abs(1/beta2 - 1/beta1) / abs(1/beta2)) > epsilon);
	return (1 / beta2 + lambda);
}
double* Doolittle(double *Y, double lambda) {
	double DD[6][502];
	double CC[5][501];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 501; j++) {
			if (i == 2) {
				CC[i][j] = C[i][j] - lambda;
			}
			else {
				CC[i][j] = C[i][j];
			}
		}
	}
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 501; j++) {
			DD[i + 1][j + 1] = CC[i][j];
		}
	}
	for (int k = 1; k <= 501; k++) {
		for (int j = k; j <= min(k + 2, 501); j++) {
			double cc = 0;
			for (int t = max(1, k - 2, j - 2); t <= k - 1; t++) {
				cc += DD[k - t + 3][t] * DD[t - j + 3][j];
			}
			DD[k - j + 3][j] -= cc;
			//cout << " DD[" << k - j + 3 << "][" << j << "] =" << DD[k - j + 3][j] << endl;
		}
		if (k != 501) {
			for (int i = k + 1; i <= min(k + 2, 501); i++) {
				double cc = 0;			
				for (int t = max(1, i - 2, k - 2); t <= k - 1; t++) {				
					cc += DD[i - t + 3][t] * DD[t - k + 3][k];			
				}
				DD[i - k + 3][k] = (DD[i - k + 3][k] - cc) / DD[3][k];
			//cout << "DD[" << i - k + 3 << "][" << k << "] =" << DD[i - k + 3][k] << endl;		
			}
		}
	}
	//回带求出X的值
	double X[501];
	double B[501];
	for (int i = 0; i < 501; i++) {
		B[i] = Y[i];
	}
	for(int i=2;i<=501;i++){
		double cc=0;
		for (int t = max(1, i - 2); t <= i - 1; t++) {
			cc += DD[i - t + 3][t] * B[t-1];
		}
		B[i-1] -= cc;
	}
	X[500] = B[500] / DD[3][501];
	for (int i = 500; i >= 1; i--) {
		double cc = 0;
		for (int t = i + 1; t <= min(i + 2, 501); t++) {
			cc += DD[i - t + 3][t] * X[t - 1];
		}
		X[i - 1] = (B[i-1] - cc) / DD[3][i];
	}
	return X;
}
double DoolittleDet() {
	double DD[6][502];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 501; j++) {
			DD[i + 1][j + 1] = C[i][j];
		}
	}
	for (int k = 1; k <= 501; k++) {
		for (int j = k; j <= min(k + 2, 501); j++) {
			double cc = 0;
			for (int t = max(1, k - 2, j - 2); t <= k - 1; t++) {
				cc += DD[k - t + 3][t] * DD[t - j + 3][j];
			}
			DD[k - j + 3][j] -= cc;
		}
		if (k != 501) {
			for (int i = k + 1; i <= min(k + 2, 501); i++) {
				double cc = 0;
				for (int t = max(1, i - 2, k - 2); t <= k - 1; t++) {
					cc += DD[i - t + 3][t] * DD[t - k + 3][k];
				}
				DD[i - k + 3][k] = (DD[i - k + 3][k] - cc) / DD[3][k];
			}
		}
	}
	double cond = 1;
	for (int i = 1; i < 502; i++) {
		cond *= DD[3][i];
	}
	return cond;
}
int max(int x, int y, int z) {
	x = x > y ? x : y;
	z = z > x ? z : x;
	return z;
}
int max(int x, int y) {
	x = x > y ? x : y;
	return x;
}
int min(int x, int y) {
	return x = x < y ? x : y;
}
