// Calc.h
// Basic Sampling and Scattering  Calculation
// Zhuoran Han
//

#ifndef CALC_H
#define CALC_H

#define pi 3.14159265358979323846

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cstdlib>

//Generate Unifrom Distribution (0,1)
double Uni_dis(); 
//Exponentail Distribution -ln(x)/Sigma   
//Sampling Walking Distance
//or Delayed Neutron (not implemented)
double Exp_dis(double Sigma); 
//Sample gamma from 0 to 2pi   
double Circ_dis();
//Sample gamma from 0 to pi   
double SemiCirc_dis();   
// G1 sampling for target nuclei velocity
//g1(x) = 4/sqrt(pi)*x^2*exp(-x^2)
double G1_dis();
// G2 sampling for target nuclei velocity
//g2(x) = 2*x^3*exp(-x^2)       
double G2_dis();
// Gmu sampling for target nuclei velocity
// u(x) =  2x-1
double Mu_dis();
// Fisson Energy Spectrum
// Watt Distribution
double Watt_dis(); 
//Sample Nu for # of fission neutron
int Nu_dis();
  

#endif