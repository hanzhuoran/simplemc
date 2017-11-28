// xsection.h
// xsections for different energy, emperical results
// Should change this part to temperature related in 
// Zhuoran Han
//
#include "xsection.h"
#include<iostream>
// Microscopic cross scetions in barn
//Calculate the cross section. c1,c2,c3 are emperical constants
double Sigma(double c1, double c2, double c3, double E)
{
	return (c1+c2/sqrt(E))*exp(c3*sqrt(E));
}

//Calculate the resonances cross section for U-238
double Sigma_res(double c1, double c2, double c3, double E)
{
	double y;
    y = 2/c3*(E-c1);
    return c2*sqrt(c1/E)/(1+y*y);
}

//Calculate the microscopic cross section for non fuel i=1-Scatter or 2-Capture or 3-total
//Hydrogen
double SigmaH(int i, double E)
{
	double s;
    switch (i) 
    {
        case 1:
            s = Sigma(2e1,3e-3,-1.2,E);
            break;
        case 2:
            s = Sigma(0,8e-5,0,E);
            break;
        case 3:
            s = Sigma(2e1,3e-3,-1.2,E)+Sigma(0,8e-5,0,E);
            break;
        default:
            break;
    }
    return s;
}

//Oxygen
double SigmaO(int i, double E)
{
	double s;
    switch (i) 
    {
        case 1:
            s = Sigma(4,1.5e-4,-6e-1,E);
            break;
        case 2:
            s = Sigma(0,0,0,E);
            break;
        case 3:
            s = Sigma(4,1.5e-4,-6e-1,E)+Sigma(0,0,0,E);
            break;
        default:
            break;
    }
    return s;
}


//Calculate the microscopic cross section for fuel. i=1-Scatter or 2-Capture or 3-Fission or 4-total
//U235
double SigmaU235(int i, double E)
{
	double s;
    switch (i) 
    {
        case 1:
            s = Sigma(1.5e1,1.5e-4,-4e-1,E);
            break;
        case 2:
            s = Sigma(4e-1,2.5e-3,-1,E);
            break;
        case 3:
            s = Sigma(8e-1,6.0e-2,0,E);
            break;
        case 4:
            s = Sigma(1.5e1,1.5e-4,-4e-1,E)
            +Sigma(4e-1,2.5e-3,-1,E)
            +Sigma(8e-1,6.0e-2,0,E);
            break;
        default:
            break;
    }
    return s;
}

//Calculate the microscopic cross section. i=1-Scatter or 2-Capture or 3-total
//U238
double SigmaU238(int i, double E)
{
	double s;
    switch (i) 
    {
        case 1:
            s = Sigma(9,1e-4,-1.6e-1,E);
            break;
        case 2:
            s = Sigma(1.8,4e-4,-1.5,E)
            +Sigma_res(6.6e-6,7e3,4e-8,E)
            +Sigma_res(2.2e-5,6e3,3e-8,E)
            +Sigma_res(3.8e-5,6.5e3,1e-7,E);
            break;
        case 3:
            s = Sigma(9,1e-4,-1.6e-1,E)
            +Sigma(1.8,4e-4,-1.5,E)
            +Sigma_res(6.6e-6,7e3,4e-8,E)
            +Sigma_res(2.2e-5,6e3,3e-8,E)
            +Sigma_res(3.8e-5,6.5e3,1e-7,E);
            break;
        default:
            break;
    }
    return s;
}




//Calculate the Macroscopic cross section in Moderator. i= 0-H, 1-O, 4-total
double Sigma_M(int i, double E)
{
	double s;
    switch (i) 
    {
        case 0:
            s = NH2O_H*SigmaH(3,E);
            break;
        case 1:
            s = NH2O_O*SigmaO(3,E);
            break;
        case 4:
            s = NH2O_H*SigmaH(3,E)+  NH2O_O*SigmaO(3,E);
            break;
        default:
            break;
    }
    return s;
}

//Calculate the Macroscopic cross section in Fuel. i= 1-O, 2-U235, 3-U238, 4-total
double Sigma_F(int i, double E, double coeff)
{
	double s;
    double f = coeff*coeff*coeff;
    switch (i)
    {
        case 1:
            s = NUO2_O/f*SigmaO(3,E);
            break;
        case 2:
            s = NUO2_U235/f*SigmaU235(4,E);
            break;
        case 3:
            s = NUO2_U238/f*SigmaU238(3,E);
            break;
        case 4:
            s = NUO2_O/f*SigmaO(3,E)+
            NUO2_U235/f*SigmaU235(4,E)+NUO2_U238/f*SigmaU238(3,E);
            break;
        default:
            break;
    }
    return s;
}

//Sample the isotope
//
int Col_iso(Neutron n, double coeff)
{
	int i;
	double sO,sH,sU235;
	double r = 0;
    double E = n.getE();
//std::cout<<"iso"<<std::endl;
	switch (n.getregion()) 
	{
//std::cout<<"region"<<n.getregion()<<std::endl;
        case 1: // Fuel
        {
//std::cout<<"in fuel "<<std::endl;
 			sO = Sigma_F(1,E,coeff)/Sigma_F(4,E,coeff);
//std::cout<<"sO"<<sO<<std::endl;
 			sU235 = (Sigma_F(1,E,coeff)+Sigma_F(2,E,coeff))/Sigma_F(4,E,coeff);
//std::cout<<"sU235"<<sU235<<std::endl;
 			r = Uni_dis();
//std::cout<<"sample"<<r<<std::endl;
 			if (r < sO) {i = 1;} //O
			else if (r < sU235) {i = 2;} //U235
 			else  {i = 3;} //U238
 			break;
        }
        case -1: // Mod
        {
//std::cout<<"in mod "<<std::endl;
 			sH = Sigma_M(0,E)/Sigma_M(4,E);
//std::cout<<"sO"<<sO<<std::endl; 			
 			r = Uni_dis();
//std::cout<<"sample"<<r<<std::endl;
 			if (r < sH) {i = 0;} // H
 			else {i = 1;} // O
            break;
        }
        default:
            break;
    }
    return i;
}  
//Sample the reaction. i is isotope. Return 1-Scattering or 2-Capture or 3-Fission
int Col_rea(Neutron n, int i)
{
	int reaction = 0;
	double E = n.getE();
	double s1,s2;
    double r = 0;
    switch (i) 
    {
        case 0: // H
        {
        	s1  = SigmaH(1,E)/SigmaH(3,E);
 			r = Uni_dis();
 			if (r < s1) {reaction = 1;} // Scattering
 			else {reaction = 2;} // Capture
            break;
        }
        case 1: // O
        {
        	s1  = SigmaO(1,E)/SigmaO(3,E);
 			r = Uni_dis();
 			if (r < s1) {reaction = 1;} // Scattering
 			else {reaction = 2;} // Capture
            break;
        }
        case 2:
        {
        	s1  = SigmaU235(1,E)/SigmaU235(4,E);
        	s2  = (SigmaU235(1,E)+SigmaU235(2,E))/SigmaU235(4,E);
 			r = Uni_dis();
 			if (r < s1) {reaction = 1;} // Scattering
 			else if (r < s2) {reaction = 2;} // Capture
 			else {reaction = 3;} //Fission 
            break;
        }
        case 3:
        {
        	s1  = SigmaU238(1,E)/SigmaU238(3,E);
 			r = Uni_dis();
 			if (r < s1) {reaction = 1;} // Scattering
 			else {reaction = 2;} // Capture
            break;
        }
        default:
            break;
    }
//std::cout<<"reaction"<<reaction<<std::endl;
    return reaction;
}