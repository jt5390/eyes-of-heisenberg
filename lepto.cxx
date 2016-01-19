#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


#include "epsilon.h"
#include "washout.h"
#include "dispersion.h"

using boost::multiprecision::cpp_dec_float_50;




int main(int argc, char * argv[])
{

    
double cutoff = 0.0;
    double M = 1e13;
bool boolMe = true;

	int c;
	while ((c = getopt (argc, argv, "WC:M:")) != -1)
   	{ 
	switch (c) 
      	{
      		case 'C':
			cutoff =strtof(optarg, NULL); 
        		break;
      		case 'M':
			M=strtof(optarg,NULL);
			break;
		case 'W':
			boolMe=true;
			break;
      		
      		case '?':
			printf("Allowed arguments:\n"
				"\t-C\tGive the lower bound cutoff for integral in epsilon_eff, Default 0.0.\n"		
				"\t-M\tInput sterile Mass in GeV, default 1e11.\n"
				"\t-W\tWashout run, default false, toggles on.\n"
				);
          	return 1;
      		default:
			printf("I don't know how you got here.\n");
        		abort ();
 	}}


	
	double baseT = 1;

if(boolMe){

   // for(double x = 1.0;x >0.01 ; x=x-0.001)
    //for(double x = 1.0;x >0.01 ; x=x-0.001)
    for(double x = 0.85;x >-0.3 ; x=x-0.01)
    {
        //double k1 =x*baseT;
        double k1=0.1*baseT;
        double k0=x*baseT;
        // double dispEquation(k0,&q);

        std::cout<<k0<<" "<<(1+sigmaA(baseT,1.5*baseT,k1,k0)*k0+sigmaB(baseT,1.5*baseT,k0,k1)-k1*(1+sigmaA(baseT,1.5*baseT,k0,k1)))<<" "<<std::endl;
       
       // double k0_09 = dispSolved(baseT, 1.5*baseT, k1);
       // double sigA = sigmaA(baseT,  M, k0_09, k1);
        //double sigB = sigmaB(baseT,  M, k0_09, k1);
//        double disp = (1+sigA)*k0_09+sigB-(1+sigA)*k1;
      //  std::cout<<sigA<<" "<<sigB<<" "<<disp<<std::endl;
       
       // std::cout<<(fabs(k1))/(baseT)<<" "<<(k0_09-fabs(k1))/(2*baseT)<<std::endl;
 

        //std::cout<<(fabs(k1))/(baseT)<<std::endl;
    }

    
    


}
else {

	//for(double z = -2; z<1; z=z+0.05){
	//	std::cout<<pow(10,z)<<" "<<epsilon_eff_ratio(M/(pow(10,z)),M,cutoff)-1.2<<std::endl;
	//}
    double ik0 = 1.2;
    baseT=1e14;
    struct dispParams2 mypar ={baseT,1.0*baseT,0.4*baseT};
    for(ik0=-0.9; ik0<0.0;ik0=ik0+0.01){
        std::cout<<pow(10,ik0)<<" wrong "<<dispEquation(pow(10,ik0)*baseT,&mypar)/baseT<<std::endl;
        
        }
    }
    
    

return 0;
}
