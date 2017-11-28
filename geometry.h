// geometry.h
// Geometry Calculation
// Zhuoran Han
//
#ifndef Geometry_H
#define Geometry_H

#include <stdio.h>
#include <math.h>
#include "neutron.h"
#include "xsection.h"

//Distance to a plane perpendicular to x at x0, given x, u
double dis_x_plane(double x, double u, double x0); 
//Distance to a plane perpendicular to y at y0, given y, v
double dis_y_plane(double y, double v, double y0); 
//Distance to a plane perpendicular to z at z0, given z, w
double dis_z_plane(double z, double w, double z0); 
//Distance to a sphere
double dis_z_sphere(double x, double y, double z, double u, double v, 
	double w, double x0, double y0, double z0, double R); 
//Distance to collison
double dis_collision(Neutron n, double factor, double coeff); 
//Minimum distance to a surface.
double dis_min(Neutron n, double x0, double y0, double z0, double R);
// Check if the neutron is out of the box, reflection boudary condition?
void boundary(Neutron &n, double x0, double y0, double z0); 
// Set the region of a neutron
void setRregion(Neutron &n, double R);      
#endif 