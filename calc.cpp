// Calc.cpp
// Basic Sampling and Scattering  Calculation
// Implementations 
// Zhuoran Han
//

#include "calc.h"

//Generate Unifrom Distribution (0,1)
double Uni_dis()
{
	//using default random number
// How to impement the random of G+1/Mod
//	srand((unsigned)time(NULL));
	return (double)rand()/(double)RAND_MAX;
}

//Exponentail Distribution -ln(x)/Sigma   
//Sampling Walking Distance
//or Delayed Neutron (not implemented)
double Exp_dis(double Sigma)
{
	return -log(Uni_dis())/Sigma;
}

//Return sampling angle in o to 2pi
double Circ_dis()
{
	return Uni_dis()*2*pi;
}
// G1 sampling for target nuclei velocity
//g1(x) = 4/sqrt(pi)*x^2*exp(-x^2)
double G1_dis()
{
	return sqrt(-log(Uni_dis())-log(Uni_dis())*pow(cos(pi/2*Uni_dis()),2));
}
// G2 sampling for target nuclei velocity
//g2(x) = 2*x^3*exp(-x^2)       
double G2_dis()
{
	return sqrt(-log(Uni_dis()*Uni_dis()));	
}
// Gmu sampling for target nuclei velocity
// u(x) =  2x-1
double Mu_dis()
{
	return 2*Uni_dis()-1;
}
// Fisson Energy Spectrum
// Watt Distribution
double Watt_dis()
{
	double a = 0.988, b = 2.249;
    double W;
    W = a*G1_dis();
    return W+a*a*b/4+(2*Uni_dis()-1)*sqrt(a*a*b*W);
}
//Sample Nu for # of fission neutron
int Nu_dis()
{
	return (int)(2.45+Uni_dis());
}