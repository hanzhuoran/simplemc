// neutron.cpp
// Neutron Class Implementation
// Zhuoran Han
//	

#include "neutron.h"

// Self Constructor
Neutron::Neutron()
{
    double r, theta, phi;
    r = Uni_dis()*D/2;
    theta = SemiCirc_dis();
    phi = Circ_dis();
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
    double mu, gamma;
    mu = Mu_dis();
    gamma = Uni_dis()*2*pi;
    Omegax = sqrt(1-mu*mu)*cos(gamma);
    Omegay = sqrt(1-mu*mu)*sin(gamma);
    Omegaz = mu;
    double newE;
    newE = Watt_dis();
    E = newE;
    T =  Troom;
    weight = 1;
    region = 1;
}
// Construct with coeff
Neutron::Neutron(double coeff)
{
    double r, theta, phi;
    r = Uni_dis()*D/2;
    theta = SemiCirc_dis();
    phi = Circ_dis();
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
    x = x*coeff;
    y = y*coeff;
    z = z*coeff;
    double mu, gamma;
    mu = Mu_dis();
    gamma = Uni_dis()*2*pi;
    Omegax = sqrt(1-mu*mu)*cos(gamma);
    Omegay = sqrt(1-mu*mu)*sin(gamma);
    Omegaz = mu;
    double newE;
    newE = Watt_dis();
    E = newE;
    T =  Troom;
    weight = 1;
    region = 1;
}
// Constructor
Neutron::Neutron(double newx,double newy,double newz,
	double newE,double newT,
	double newOmegax,double newOmegay,double newOmegaz)
{
	x = newx; y = newy; z = newz;
    E = newE; T = newT;
    Omegax = newOmegax; Omegay = newOmegay; Omegaz = newOmegaz;
    weight = 1; // new neutron has weight as 1
    region = 1; // new neutron is born in fuel region
}

//get basic info
double Neutron::getX() {return x;}
double Neutron::getY() {return y;}
double Neutron::getZ() {return z;}
double Neutron::getE() {return E;}
double Neutron::getT() {return T;}
double Neutron::getOmegax() {return Omegax;}
double Neutron::getOmegay() {return Omegay;}
double Neutron::getOmegaz() {return Omegaz;}

// get weight, return 0(not exist) or 1(exist)
int Neutron::getweight() {return weight;} 
// get region, return 1 or -1
int Neutron::getregion() {return region;} 

// set region. 1 for Fuel, -1 for mod
void Neutron::setregion(int r) {region = r;}
// set position
void Neutron::setposition(double newx,double newy,double newz)
{
	x = newx; 
	y = newy; 
	z = newz;
} 
// set velocity
void Neutron::setvelocity(double newOmegax,double newOmegay,double newOmegaz)
{
	Omegax = newOmegax;
	Omegay = newOmegay; 
	Omegaz = newOmegaz;
}
//capture, weight = 0, vanish
void Neutron::capture()
{
	weight = 0;
}
//i indicates nuclei species: 0-H, 1-O, 2-U235,3-U238
void Neutron::scattering(int i)
{
    double beta, x, y, omega1, V_tilde, mu_tilde = 1, v, 
    	Omegax_T, Omegay_T, Omegaz_T, gamma_T;
    v = sqrt(2*E/m_n);
    
    beta = sqrt(m_n*A[i]/(2*k_B*T));
    y = beta*v;
    omega1 = sqrt(pi)*y/(2+sqrt(pi)*y);
    V_tilde = v;
    while (Uni_dis() > sqrt(v*v+V_tilde*V_tilde-2*v*V_tilde*mu_tilde)/(v+V_tilde)) {
        if (Uni_dis() < omega1) {
            x = G1_dis();
        } else {
            x = G2_dis();
        }
        V_tilde = x/beta;
        mu_tilde = Mu_dis();
    }
    gamma_T = Circ_dis();
    Omegax_T = mu_tilde*Omegax+(Omegax*Omegaz*cos(gamma_T)-Omegay*sin(gamma_T))*sqrt((1-mu_tilde*mu_tilde)/(1-Omegaz*Omegaz));
    Omegay_T = mu_tilde*Omegay+(Omegay*Omegaz*cos(gamma_T)-Omegax*sin(gamma_T))*sqrt((1-mu_tilde*mu_tilde)/(1-Omegaz*Omegaz));
    Omegaz_T = mu_tilde*Omegaz-cos(gamma_T)*sqrt((1-mu_tilde*mu_tilde)*(1-Omegaz*Omegaz));
    
    
    double u_x, u_y, u_z;
    u_x = (v*Omegax + A[i]*V_tilde*Omegax_T)/(1+A[i]);
    u_y = (v*Omegay + A[i]*V_tilde*Omegay_T)/(1+A[i]);
    u_z = (v*Omegaz + A[i]*V_tilde*Omegaz_T)/(1+A[i]);
    
    double vc_x, vc_y, vc_z, vc;
    double Omegax_c, Omegay_c, Omegaz_c, mu_c, gamma_c; 
    double Omegax_cprime, Omegay_cprime, Omegaz_cprime;
    vc_x = v*Omegax-u_x;
    vc_y = v*Omegay-u_y;
    vc_z = v*Omegaz-u_z;
    
    vc = sqrt(vc_x*vc_x+vc_y*vc_y+vc_z*vc_z);
    Omegax_c = vc_x/vc;
    Omegay_c = vc_y/vc;
    Omegaz_c = vc_z/vc;
    
    mu_c = Mu_dis();
    gamma_c = Circ_dis();
    
    Omegax_cprime = mu_c*Omegax_c+(Omegax_c*Omegaz_c*cos(gamma_c)-Omegay_c*sin(gamma_c))*sqrt((1-mu_c*mu_c)/(1-Omegaz_c*Omegaz_c));
    Omegay_cprime = mu_c*Omegay_c+(Omegay_c*Omegaz_c*cos(gamma_c)-Omegax_c*sin(gamma_c))*sqrt((1-mu_c*mu_c)/(1-Omegaz_c*Omegaz_c));
    Omegaz_cprime = mu_c*Omegaz_c-cos(gamma_c)*sqrt((1-mu_c*mu_c)*(1-Omegaz_c*Omegaz_c));
    
    double v_labprime_x, v_labprime_y, v_labprime_z, v_labprime;
    v_labprime_x = vc*Omegax_cprime;
    v_labprime_y = vc*Omegay_cprime;
    v_labprime_z = vc*Omegaz_cprime;
    
    v_labprime = sqrt(v_labprime_x*v_labprime_x+v_labprime_y*v_labprime_y+v_labprime_z*v_labprime_z);
    
    E = m_n*v_labprime*v_labprime/2;
    Omegax = v_labprime_x/v_labprime;
    Omegay = v_labprime_y/v_labprime;
    Omegaz = v_labprime_z/v_labprime;
}

//fission, old neutron vanishes, weight = 0//?implementing generation of neutrons here?
void Neutron::fission()
{
	weight = 0;
}
// Leakage from z direction
void Neutron::leak()
{
	weight = 0;
}
