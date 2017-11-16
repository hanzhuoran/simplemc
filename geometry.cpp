// geometry.cpp
// Geometry Calculation Implementaion
// Zhuoran Han
//

#include "geometry.h"

//Distance to a plane perpendicular to x at x0, given x, u
double dis_x_plane(double x, double u, double x0)
{
	return (x0-x)/u;
}

//Distance to a plane perpendicular to y at y0, given y, v
double dis_y_plane(double y, double v, double y0)
{
	return (y0-y)/v;
}

//Distance to a plane perpendicular to z at z0, given z, w
double dis_z_plane(double z, double w, double z0)
{
	return (z0-z)/w;
}

//Distance to a cylinder parallel to z-axis
double dis_z_cylinder(double x, double y, double u, double v, 
	double x0, double y0, double R)
{
	double x1, y1, a, k, c, d;
    x1 = x-x0;
    y1 = y-y0;
    a = u*u+v*v;
    k = x1*u+y1*v;
    c = x1*x1+y1*y1-R*R;
    if (a == 0 || k*k-a*c < 0) 
    {
        d = 1000;      
        //No intersections with the cylinder. set d to a very large value
    } 
    else if (c < 0) 
    {
        d = (-k+sqrt(k*k-a*c))/a;
    } 
    // c >=0
    else if ((-k-sqrt(k*k-a*c)) > 0) 
    {
        d = (-k-sqrt(k*k-a*c))/a;
    }
    else 
    {
        d = 1000;       
        //default is a very large number
    }       
    return d;
}

//Distance to collison
double dis_collision(Neutron n)
{
	double d;
	double f = factor*factor*factor;
    switch (n.getregion()) 
    {
        case 1:
        {
            d = Exp_dis(Sigma_F(4, n.getE())/f);
            break;
        }
        case -1:
        {
        	//we don't expand water
            d = Exp_dis(Sigma_M(4, n.getE())/f);
            break;
        }
        default:
            break;
    }
    return d;
}

//Minimum distance to a surface.
double dis_min(Neutron n, double x0, double y0, double z0, double R)
{
	double d[7];
	double m = 10000; //very large number
    d[0] = dis_x_plane(n.getX(), n.getOmegax(), x0);
    d[1] = dis_x_plane(n.getX(), n.getOmegax(), -x0);
    d[2] = dis_y_plane(n.getY(), n.getOmegay(), y0);
    d[3] = dis_y_plane(n.getY(), n.getOmegay(), -y0);
    d[4] = dis_z_plane(n.getZ(), n.getOmegaz(), z0);
    d[5] = dis_z_plane(n.getZ(), n.getOmegaz(), -z0);
    d[6] = dis_z_cylinder(n.getX(), n.getY(), n.getOmegax(), n.getOmegay(), 0, 0, R);
    for (int i = 0; i < 7; i++) 
    {
        if (d[i] > 0 && d[i] < m) //if negative, neutron goes away
        {
            m = d[i];
        }
    }
    return m;
}

// Check if the neutron is out of the box.
// Reflective boudary condition on X and Y. Leakage on Z.
void boundary(Neutron &n, double x0, double y0, double z0)
{
	if (n.getX() >= x0)
	{
        n.setposition(2*x0-n.getX(), n.getY(), n.getZ());
        n.setvelocity(-n.getOmegax(), n.getOmegay(), n.getOmegaz());
    }
    if (n.getX() <= -x0) 
    {
        n.setposition(-2*x0-n.getX(), n.getY(), n.getZ());
        n.setvelocity(-n.getOmegax(), n.getOmegay(), n.getOmegaz());
    }
    if (n.getY() >= y0) 
    {
        n.setposition(n.getX(), 2*y0-n.getY(), n.getZ());
        n.setvelocity(n.getOmegax(), -n.getOmegay(), n.getOmegaz());
    }
    if (n.getY() <= -y0) 
    {
        n.setposition(n.getX(), -2*y0-n.getY(), n.getZ());
        n.setvelocity(n.getOmegax(), -n.getOmegay(), n.getOmegaz());
    }
    if (n.getZ() >= z0) 
    {
    	n.leak();
//std::cout<<"Leakage"<<std::endl;
    }
    if (n.getZ() <= -z0)
    {
    	n.leak();
//std::cout<<"Leakage"<<std::endl;
    }
}

// Set the region of a neutron based on fuel radius
void setRregion(Neutron &n, double R)
{
	if (n.getX()*n.getX()+n.getY()*n.getY() < R*R) {n.setregion(1);} 
    else {n.setregion(-1);}
}
