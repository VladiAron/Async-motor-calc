#include <stdio.h>
#include <math.h>

//page 425 Ток I0
//page 426 theory
//page 472 practice

#define HEADER "\"S\", \"P2\", \"P1\", \"КПД\", \"COS\", \"I1\", \"M\", \"NR\",\n"
#define OUTPUT_FORMAT "%'.3f, %'.3f, %'.3f, %'.3f, %'.3f, %'.3f, %'.3f, %'.1f,\n"


double s[30];
double P2[30];
double P1[30];
double KPD[30];
double cosfi[30];
double I1[30];
double M[30];
double NR[30];


double U1 = 220, p = 2, f = 50,
		x1 = 29.17, r1 = 53.74, _x2 = 35.58,
		_r2 = 50, x12 = 378.17, r12 = 17.7,
		hnr = 0.0083, myr = 1, ro = 0.044/1000000,
		W1 = 720, kob = 0.97, tau = 0.0432, kmy = 1.44;




int main(int argc, char * argv[]){

	double 	I0_a = 0.52, I0_p = 7.91;
	double prir = 0.005;
	double c1, c1_a, c1_p,_a, a, _b, b;

	double n = 60 * f / p;
//	c1_a = (r12 * (r1 + r12) + x1 * (x1 + x12)) / (pow(r12, 2) + pow(x12, 2));
//	c1_p = (x1 * r12 - r1 * x12) / (pow(r12, 2) + pow(x12, 2));
//
//	c1 = sqrt(pow(c1_a, 2) + pow(c1_p, 2));
//
//	_a = pow(c1_a, 2) - pow(c1_p, 2);
//	_b = 2 * c1_a * c1_p;
//	a = c1_a * r1 - c1_p * x1 - _b * _x2;
//	b = c1_a * x1 + c1_p * r1 + _a * _x2;

	int i = 0;

	s[0] = 0.005;

	double _I2;
	double Z12;
	double I2n;
	double _Z2;
	double E1;
	double Imy = 1;
	double Z1;
	double Up = 0;
	double x12n = 1;


	do{

		NR[i] = n * (1 - s[i]);
		//M[i] = (3 * pow(U1, 2) * p * _r2)/(2 * M_PI * f * s[i] *(pow(r1 + _r2/s[i], 2) + pow(x1 + _x2, 2)));

//		double res1 = _a * _r2 / s[i];
//		double res2 = _b * _r2 / s[i];
//		double R = a + res1;
//		double X = b + res2;
//		double Z = sqrt(pow(R, 2) + pow(X ,2));
//		double __I2 = U1 / Z;
//		double _cosfi2 = R / Z;
//		double _sinfi2 = X / Z;
//		double I1_a = I0_a + __I2 * _cosfi2;
//		double I1_p = I0_p + __I2 * _sinfi2;
//		I1[i] = sqrt(pow(I1_a, 2) + pow(I1_p, 2));
//		_I2 = c1 * __I2;
//		double I1_a = I1[i] * cosfi[i];
//		P1[i] = 3 * U1 * I1_a / 1000;
		//cosfi[i] = I1_a / I1[i];

		_I2 = 0;
		Up = 0;
		while (Up < U1){
			_I2 += 0.001;
			_Z2 = sqrt(pow(_x2, 2) + pow(_r2/s[i], 2));
			E1 = _I2 * _Z2;
			//x12n = x12 * kmy - (kmy - 1)*(x12 - x1) * E1 / U1;
			Z12 = sqrt(pow(x12, 2) + pow(r12, 2));
			Imy = E1/Z12;
			I1[i] = Imy + _I2;
			Z1 = sqrt(pow(x1, 2) + pow(r1, 2));
			Up = E1 + I1[i] * Z1;
		}
		double Pm = 3 * pow(Imy,2) * r12;
		double Pe1 = 3 * pow(I1[i], 2) * r1;
		double Pem = 3 * pow(_I2, 2) * _r2/s[i];
		P1[i] = Pem + Pm + Pe1;
		double Pe2 = 3 * pow(_I2, 2) * _r2;
		double Pmex = 3 * pow(_I2,2) * _r2 * (1 - s[i])/s[i];//Pem - Pe2;
		M[i] = Pmex * p / (2 * M_PI * f * (1 - s[i]));
		double Pdob = 0.005 * P1[i];
		P2[i] = Pmex - Pdob;
		double EP = Pe2 + Pe1 + Pdob;
		KPD[i] = (P2[i] / P1[i]);
		double cosff =  (3 * pow(I1[i], 2) * Z1);
		i++;
		s[i] = s[i-1] + prir;
		if(s[i] >= 0.199){
			prir = 0.1;
		}else if(s[i] >= 0.019){
			prir = 0.02;
		}

	}while (s[i] <= 1);

	printf(HEADER);
	for(int j = 0; j < i; j++){
		printf(OUTPUT_FORMAT, s[j], P2[j], P1[j], KPD[j], KPD[j], I1[j], M[j], NR[j]);
	}



	return 0;

}




