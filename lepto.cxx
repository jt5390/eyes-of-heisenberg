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

#include "epsilon.h"
#include "washout.h"
#include "dispersion.h"

int main(int argc, char * argv[])
{

double cutoff = 0.0;
double M = 1e11;
bool boolMe = false;

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


	std::cout<<"#Cutoff: "<<cutoff<<std::endl;


if(boolMe){
	
	/*	double T = 1e12;
	struct washoutParams params = {T,M};
	for(double x = 0; x<20; x=x+0.2)
	{
		std::cout<<x<<" "<<washout_integrand(x,&params)<<std::endl;
	}
	for(double z = -2; z<1; z=z+0.05){
		std::cout<<pow(10,z)<<" "<<washout(1,M/(pow(10,z)),M)<<std::endl;
	}*/

	struct dispParams params = {1e11,0.7e11,2e12,3e12};
/*	for(double x = 10; x<14; x=x+0.0101)
	{
//		std::cout<<pow(10,x)<<" "<<TRKS_integrand(pow(10,x),&params)<<TRUS_integrand(pow(10,x),&params)<<std::endl;
		std::cout<<pow(10,x)<<" "<<dispEquation(pow(10,x),&params)<<std::endl;
	}

	std::cout<<"# integral: "<<TRKS(1e12,1e11,2e11,3e11)<<"  "<<TRUS(1e12,1e11,2e11,3e11)<<std::endl;
	std::cout<<"# a and b test: "<<sigmaA(params)<<"  "<<sigmaB(params)<<std::endl;
	std::cout<<"# disp : "<<dispSolved(1.2e12,params)<<std::endl;*/

	
	struct dispParams params2={1e11,0.9e11,2e12,3e12};
	struct dispParams params3={1e11,1e11,2e12,3e12};
	struct dispParams params4={1e11,1.2e11,2e12,3e12};
	struct dispParams params5={1e11,1.5e11,2e12,3e12};
	struct dispParams params6={1e11,2e11,2e12,3e12};

	for(double x = 0.05*params.T; x<params.T; x=x+0.01*params.T)
	{
		double k = dispSolved(x,params);
		double kB= dispSolved(x,params2);
		double kC = dispSolved(x,params3);
		double kD = dispSolved(x,params4);
		double ke = dispSolved(x,params5);
		double kf = dispSolved(x,params6);
		double kz = x;
//		std::cout<<pow(10,x)<<" "<<TRKS_integrand(pow(10,x),&params)<<TRUS_integrand(pow(10,x),&params)<<std::endl;
		std::cout<<kz/params.T<<" "<<(k-kz)/params.T<<"  "<<(kB-kz)/params.T<<"  "<<(kC-kz)/params.T<<" "<<(kD-kz)/params.T<<" "<<(ke-kz)/params.T<<"  "<<(kf-kz)/params.T<<std::endl;
	}


} else {
	/*double T = 1e12;
	struct epsilonParams params = {T,M};
	for(double kk = 4; kk<13.5; kk=kk+0.1)
	{
		std::cout<<pow(10,kk)<<" "<<epsilon_numerator_integrand (pow(10,kk),&params)<<" "<<epsilon_denominator_integrand(pow(10,kk),&params)<<std::endl;
	}*/

	for(double z = -2; z<1; z=z+0.05){
		std::cout<<pow(10,z)<<" "<<epsilon_eff_ratio(M/(pow(10,z)),M,cutoff)-1.0<<std::endl;
	}
}
return 0;
}
