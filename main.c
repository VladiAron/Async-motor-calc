#include <stdio.h>
#include <math.h>

//page 426 theory
//page 472 practice

double s[30];
double P2[30];
double P1[30];
double KPD[30];
double cosfi[30];
double I1[30];
double M[30];
double NR[30];



int main(int argc, char * argv[]){

	double P2_nom = 15, U1 = 220, p = 4,
			I0_a = 0.52, I0_p = 7.91, I_m = 7.91,
			P_ct_mex = 0.49, r1 = 0.355, _r2 = 0.186,
			c1 = 1.025, _a = 1.051, a = 0.364, _b = 0, b = 1.65, P_esh = 0;


	int i = 0;

	s[0] = 0.005;

	do{
		s[i] = 0.005 + 0.005*i;
		double res1 = _a * _r2 / s[i];
		double res2 = _b * _r2 / s[i];
		double R = a + res1;
		double X = b + res2;
		double Z = sqrt(pow(R, 2) + pow(X ,2));
		double __I2 = U1 / Z;
		double _cosfi2 = R / Z;
		double _sinfi2 = X / Z;
		double I1_a = I0_a + __I2 * _cosfi2;
		double I1_p = I0_p + __I2 * _sinfi2;
		I1[i] = sqrt(pow(I1_a, 2) + pow(I1_p, 2));
		double _I2 = c1 * __I2;
		P1[i] = 3 * U1 * I1_a / 1000;
		double Pe1 = 3 * pow(I1[i], 2) * r1 / 1000;
		double Pe2 = 3 * pow(_I2, 2) * _r2 / 1000;
		double Pdob = 0.005 * P1[i];
		double EP = P_ct_mex + Pe2 + Pe1 + Pdob + P_esh;
		P2[i] = P1[i] - EP;
		KPD[i] = 1 - (EP / P1[i]);
		cosfi[i] = I1_a / I1[i];
		i++;

	}while (s[i-1] < 0.030);


	return 0;

}
