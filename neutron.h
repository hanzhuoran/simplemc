// neutron.h
// Neutron Class
// Zhuoran Han
//

#ifndef Neutron__H
#define Neutron__H

#define Troom 298.15
#define D 3
#define L 50

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cstdlib>
#include "calc.h"

//Atom mass for 0-H, 1-O, 2-U235, 3-U238
const double A[4] = {1.00794, 15.9994, 235.0439299, 238.05078826};      
//Atom mass of neutron
const double m_n = 1.00866491600; 
//Boltzmann constant in unit (MeV/K)
const double k_B = 8.6173324E-11;
//Tempreture in
const double Tin = 293.0;

class Neutron
{
public:
	//constuctor
	Neutron();
	//new neutron
	Neutron(double newx,double newy,double newz,
		double newE,double newT,
		double newOmegax,double newOmegay,double newOmegaz);
	//get basic info
	double getX();
	double getY();
	double getZ();
	double getE();
	double getT();
	double getOmegax();
	double getOmegay();
	double getOmegaz();
	int getweight(); // get weight, return 0(not exist) or 1(exist)
	int getregion(); // get region, return 1 or -1
	//set info
	void setregion(int i);// set region. 1 for Fuel, -1 for mod
//should have a set T for temperature deformation
	void setposition(double newx,double newy,double newz); // set position
	void setvelocity(double newOmegax,double newOmegay,double newOmegaz); // set velocity
	
	void capture(); //capture, weight = 0, vanish
	void scattering(int i);//i indicates nuclei species: 0-H, 1-O, 2-U235,3-U238
	void fission(); //fission, old neutron vanishes, generate new neutron
	void leak(); // Leakage from z direction
private:
	double x,y,z;   //Position
    double E;       //Energy
    double T;		//Temperature
    double Omegax,Omegay,Omegaz;        //Direction
    int region;     //Region 1-fuel -1-moderator
    int weight;     //weight is 0 when it is absorbed, otherwise it keeps 1

};



#endif 
