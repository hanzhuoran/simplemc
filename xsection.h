// xsection.h
// xsections for different energy, emperical results
// Should change this part to temperature related in 
// Zhuoran Han
//
#ifndef Xsection_H
#define Xsection_H

// 4% of U-235 in UO2
// 1/(cm*barn)
#define NH2O_H 6.6911e-2 
#define NH2O_O 3.3455e-2   
#define NUO2_O 4.7284e-2   
#define NUO2_U235 9.4567e-4  
#define NUO2_U238 2.2696e-2  

#include<stdio.h>
#include<iostream>
#include<math.h>
#include "neutron.h"
#include "calc.h"

//Calculate the cross section. c1,c2,c3 are emperical constants
double Sigma(double c1, double c2, double c3, double E); 
//Calculate the resonances cross section for U-238
double Sigma_res(double c1, double c2, double c3, double E); 
//Calculate the microscopic cross section for non fuel i=1-Scatter or 2-Capture or 3-total
//Hydrogen
double SigmaH(int i, double E); 
//Oxygen
double SigmaO(int i, double E); 
//Calculate the microscopic cross section for fuel. i=1-Scatter or 2-Capture or 3-Fission or 4-total
//U235
double SigmaU235(int i, double E);        
//U238
double SigmaU238(int i, double E); 
//Calculate the Macroscopic cross section in Fuel. i=0-H, 1-O, 2-U235, 3-U238, 4-total
double Sigma_F(int i, double E);        
//Calculate the Macroscopic cross section in Moderator. i=0-H, 1-O, 2-U235, 3-U238, 4-total
double Sigma_M(int i, double E);       
//Sample the isotope
int Col_iso(Neutron n);   
//Sample the reaction i-isotope return 1-Scatter or 2-Capture or 3-Fission
int Col_rea(Neutron n, int i);      
#endif 