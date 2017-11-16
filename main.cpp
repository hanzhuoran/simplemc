// main.cpp
// Simple MC codes to evalute Virtual Density
// Zhuoran Han

#include <math.h>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "neutron.h"
#include "calc.h"
#include "xsection.h"
#include "geometry.h"

const int source_number=100000;
const int iteration=40;
const double delta_d = 1E-15;
//const double Trm = 298.15;

using namespace std;

int main ()
{
	srand((unsigned)time(0));
	vector <Neutron> sourceBank(source_number); //source
    vector <Neutron> fissionBank;
//cout <<"X" << fissionBank[0].getX() << endl;
//cout <<"Y" << fissionBank[0].getY() << endl;
//cout <<"Z" << fissionBank[0].getZ() << endl;
    double k[iteration]= {0};
    double pitch = 4.0; 
    double R = D/2;
    double x0, y0, z0;
    x0 = y0 = pitch/2;
    z0 = L/2;
    R = R*coeff;
    x0 = x0*coeff;
    y0 = y0*coeff;
    z0 = z0*coeff;
    double dis_surface, dis_coll;
    dis_surface = dis_coll = 0;
    int isotope;
    int reaction;
    int nu;
    double average_k = 0;
    //End of Declaration
    for (int i = 0; i < iteration; i++)
    {
//cout << "new "<<sourceBank[0].getweight()<<endl;
//cout << "size"<<sourceBank.size()<<endl;
//cout<<endl;
//cout << "iter"<<i<<endl;

    	for (int j = 0; j < sourceBank.size(); j++)
    	{
//cout << "start "<<sourceBank[j].getweight()<<endl;
    		while(sourceBank[j].getweight()!=0)
    		{

				dis_surface = dis_min(sourceBank[j], x0, y0, z0, R);
                dis_coll = dis_collision(sourceBank[j]);
//cout <<"X" << sourceBank[j].getX() << endl;
//cout <<"Y" << sourceBank[j].getY() << endl;
//cout <<"Z" << sourceBank[j].getZ() << endl;
//cout << dis_coll <<"dis" <<dis_surface <<endl;                
                if(dis_coll < dis_surface)
                {
// cout<<"Collison"<<endl;
					sourceBank[j].setposition((sourceBank[j].getX()+sourceBank[j].getOmegax()*(dis_coll+delta_d)),
						(sourceBank[j].getY()+sourceBank[j].getOmegay()*(dis_coll+delta_d)),
						(sourceBank[j].getZ()+sourceBank[j].getOmegaz()*(dis_coll+delta_d)));
                    //Sampling isotope and reaction
                    isotope = Col_iso(sourceBank[j]); 
 					reaction = Col_rea(sourceBank[j], isotope);
//cout<< isotope <<" and " <<reaction << endl;
 					switch (reaction)
 					{     
 						case 1: //Scattering
 						{
 							sourceBank[j].scattering(isotope);

 							break;
 						}
 						case 2: // Capture
 						{
 							sourceBank[j].capture();
 							break;
 						}
 						case 3: // Fission
 						{

 							sourceBank[j].fission();
 							nu = Nu_dis();
//cout << "Fission"<<nu<<endl;
 							for (int l = 0; l < nu; l++)
 							{
                                double newmu, newgamma;
                                newmu = Mu_dis();
                                newgamma = Circ_dis();
                                Neutron f(sourceBank[j].getX(),sourceBank[j].getY(),sourceBank[j].getZ(),
                                	Watt_dis(), Troom,
                                	sqrt(1-newmu*newmu)*cos(newgamma),sqrt(1-newmu*newmu)*sin(newgamma),newmu);
                                fissionBank.push_back(f);
 							}
 							break;
 						}
 						default:
 							break;
                    }
                }
                else
                {

 					sourceBank[j].setposition(sourceBank[j].getX()+
 						sourceBank[j].getOmegax()*(dis_surface+delta_d), 
 						sourceBank[j].getY()+sourceBank[j].getOmegay()*(dis_surface+delta_d),
 						sourceBank[j].getZ()+sourceBank[j].getOmegaz()*(dis_surface+delta_d));
 					boundary(sourceBank[j], x0, y0, z0);
                    setRregion(sourceBank[j], R);
//cout<<"Stop"<<endl;
//cout << "getweight "<<sourceBank[0].getweight()<<endl; 
                }
    		}
    	}

    	//End of All History
		// should we use fissionBank_size or still fixed source_number
		long fissionBank_size;
        fissionBank_size = fissionBank.size();
//cout<<"HERE "<< fissionBank_size <<endl;
        if (fissionBank_size >= 1)
        {
        	k[i] = (double)(fissionBank.size())/(double)(sourceBank.size());
        	sourceBank.clear();
        	for (int p = 0; p < source_number; p++)
        	{
            	int index;
            	index = floor(Uni_dis()*fissionBank_size);
				sourceBank.push_back(fissionBank[index]);
        	}
        	fissionBank.clear();
        }
        else
        {
//cout << "lol" <<endl;
        	sourceBank.clear();
        	for (int p = 0; p < source_number; p++)
        	{
//cout<<"bla"<<endl;
        		Neutron n;
//cout << "new dot"<<n.getweight()<<endl;
        		sourceBank.push_back(n);
        	}

        }
 //cout << "new dot"<<sourceBank[0].getweight()<<endl;
    }
    // How to calculate k-eff? Average all batches?
    for (int s = 0; s < iteration; s++) 
    {
    	average_k = average_k + k[s]/iteration;
	}
	cout<<"keff = "<<average_k<<endl;
	return 0;
}